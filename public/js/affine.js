/*
 * 仿射拟合（浏览器 / Node 通用模块）
 *
 * 给定 n ≥ 3 对控制点 (x_i, y_i) → (u_i, v_i)，求解六参数仿射变换：
 *     u = a·x + b·y + c
 *     v = d·x + e·y + f
 *
 * 提供两种模式：
 *  1) fitAffine 经典最小二乘：使全部控制点残差平方和最小。
 *  2) fitAffineRobust 稳健模式：按固定顺序枚举有上限的非共线三点候选
 *     （三点恰好确定一个仿射模型），以「内点数量优先、内点平方误差次之」
 *     选出最佳候选，再对全部内点重新做最小二乘拟合，并迭代重分类至稳定。
 *     内点判定基于参考图像素容差 tolerancePx；外点不参与最终 RMSE。
 *     枚举上限保证点对增多时计算量有界，拖动不卡顿。
 *
 * 数值实现：最小二乘先对源点与目标点分别中心化，正规方程退化为两个
 * 2×2 线性方程组（平移量由质心回代），避免大坐标下的数值病态。
 *
 * 近共线判定：点集协方差矩阵的最小/最大特征值之比小于
 * COLLINEAR_RATIO 时视为近共线（含全部重合），拒绝拟合。
 *
 * 逐点留出验证（LOOCV）：每个参与拟合的控制点被单独移除，用其余参与点
 * 重新做最小二乘仿射拟合，再预测留出点的参考图位置并计算误差；稳健模式
 * 仅以当前内点为参与集，外点标记为「不参与」。参与点不足 4 个，或某次
 * 留出后的点集近共线时，该点标记为「无法验证」，绝不用原拟合残差充数。
 * 汇总（RMSE/最大误差/数量）仅统计可验证点，零个时为 null（不可用）。
 */
(function (root, factory) {
  const api = factory();
  if (typeof module !== 'undefined' && module.exports) module.exports = api;
  else root.Affine = api;
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  const MIN_POINTS = 3;
  // 逐点留出验证要求至少 4 个参与点：留出 1 个后仍需 ≥3 个非共线点才能重拟合
  const MIN_LOO_POINTS = 4;
  // 近共线阈值：λmin/λmax < 1e-4，即点集垂直方向展宽不足主方向的 1%
  const COLLINEAR_RATIO = 1e-4;
  // 稳健模式：枚举三点候选（及扫描三元组）的上限，保证计算量有界
  const MAX_CANDIDATES = 256;
  // 稳健模式：内点重拟合后重新分类的最大迭代次数（实测通常 1~2 次即稳定）
  const MAX_REFINE_ITERS = 16;
  // 稳健模式默认参考图像素容差
  const DEFAULT_TOLERANCE_PX = 2;

  class FitError extends Error {
    constructor(code, message) {
      super(message);
      this.name = 'FitError';
      this.code = code;
    }
  }

  function centroid(points) {
    let sx = 0, sy = 0;
    for (const p of points) { sx += p.x; sy += p.y; }
    return { x: sx / points.length, y: sy / points.length };
  }

  /**
   * 评估点集共线程度。
   * @returns {{degenerate: boolean, ratio: number}} ratio 为 λmin/λmax
   */
  function collinearity(points) {
    const c = centroid(points);
    let sxx = 0, sxy = 0, syy = 0;
    for (const p of points) {
      const dx = p.x - c.x, dy = p.y - c.y;
      sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
    }
    const trace = sxx + syy;
    if (trace <= 0) return { degenerate: true, ratio: 0 }; // 所有点重合
    let det = sxx * syy - sxy * sxy;
    if (det < 0) det = 0; // 浮点误差兜底
    const disc = Math.sqrt(Math.max(trace * trace - 4 * det, 0));
    const lMax = (trace + disc) / 2;
    const lMin = (trace - disc) / 2;
    const ratio = lMax > 0 ? lMin / lMax : 0;
    return { degenerate: ratio < COLLINEAR_RATIO, ratio };
  }

  function validate(src, dst) {
    if (!Array.isArray(src) || !Array.isArray(dst) || src.length !== dst.length) {
      throw new FitError('BAD_INPUT', '源点与目标点必须为等长数组');
    }
    if (src.length < MIN_POINTS) {
      throw new FitError('TOO_FEW_POINTS', `至少需要 ${MIN_POINTS} 对控制点`);
    }
    for (const p of src.concat(dst)) {
      if (!p || !Number.isFinite(p.x) || !Number.isFinite(p.y)) {
        throw new FitError('BAD_INPUT', '控制点坐标必须为有限数字');
      }
    }
  }

  /**
   * 最小二乘求解（不做校验）：调用方需保证至少 3 对点且两侧均非共线。
   */
  function solveAffine(src, dst) {
    const n = src.length;
    const cs = centroid(src), cd = centroid(dst);
    let sxx = 0, sxy = 0, syy = 0, sxu = 0, syu = 0, sxv = 0, syv = 0;
    for (let i = 0; i < n; i++) {
      const dx = src[i].x - cs.x, dy = src[i].y - cs.y;
      const du = dst[i].x - cd.x, dv = dst[i].y - cd.y;
      sxx += dx * dx; sxy += dx * dy; syy += dy * dy;
      sxu += dx * du; syu += dy * du;
      sxv += dx * dv; syv += dy * dv;
    }
    const D = sxx * syy - sxy * sxy; // 非共线检查已保证 D > 0

    const a = (sxu * syy - syu * sxy) / D;
    const b = (syu * sxx - sxu * sxy) / D;
    const d = (sxv * syy - syv * sxy) / D;
    const e = (syv * sxx - sxv * sxy) / D;
    const c = cd.x - a * cs.x - b * cs.y;
    const f = cd.y - d * cs.x - e * cs.y;
    return { a, b, c, d, e, f };
  }

  function residualsOf(matrix, src, dst) {
    return src.map((p, i) => {
      const u = matrix.a * p.x + matrix.b * p.y + matrix.c;
      const v = matrix.d * p.x + matrix.e * p.y + matrix.f;
      return Math.hypot(u - dst[i].x, v - dst[i].y);
    });
  }

  /**
   * 仿射最小二乘拟合（经典模式，使用全部控制点）。
   * @param {Array<{x:number,y:number}>} src 源点（扫描图）
   * @param {Array<{x:number,y:number}>} dst 目标点（参考图）
   * @returns {{matrix:{a,b,c,d,e,f}, residuals:number[], rmse:number, maxIndex:number, pointCount:number, inliers:boolean[], leaveOneOut:object}}
   * @throws {FitError} code ∈ TOO_FEW_POINTS | COLLINEAR_SRC | COLLINEAR_DST | BAD_INPUT
   */
  function fitAffine(src, dst) {
    validate(src, dst);
    const n = src.length;

    if (collinearity(src).degenerate) {
      throw new FitError('COLLINEAR_SRC', '扫描图控制点过于接近共线（或重合），无法求解仿射变换，请调整点位分布');
    }
    if (collinearity(dst).degenerate) {
      throw new FitError('COLLINEAR_DST', '参考图控制点过于接近共线（或重合），无法求解仿射变换，请调整点位分布');
    }

    const matrix = solveAffine(src, dst);
    const residuals = residualsOf(matrix, src, dst);
    const rmse = Math.sqrt(residuals.reduce((s, r) => s + r * r, 0) / n);
    let maxIndex = 0;
    residuals.forEach((r, i) => { if (r > residuals[maxIndex]) maxIndex = i; });

    // 经典模式下全部控制点均参与拟合与留出验证
    const leaveOneOut = leaveOneOutValidate(src, dst, new Array(n).fill(true));

    return { matrix, residuals, rmse, maxIndex, pointCount: n, inliers: new Array(n).fill(true), leaveOneOut };
  }

  /**
   * 稳健仿射拟合（RANSAC-lite）。
   *
   * 按 (i,j,k) 字典序枚举至多 maxCandidates 个三元组，仅评估两侧均非共线
   * 的候选；以三点精确解作为候选模型，统计参考图像素容差内的内点：
   *   内点数更多者胜；内点数相同则内点平方误差和更小者胜（首个最优点保留）。
   * 选定模型后，对内点重新做最小二乘拟合，并按新模型重新分类，迭代至
   * 内点集合稳定。最终残差对所有点给出，RMSE 仅统计内点。
   *
   * @param {Array<{x:number,y:number}>} src 源点（扫描图）
   * @param {Array<{x:number,y:number}>} dst 目标点（参考图）
   * @param {{tolerancePx?:number, maxCandidates?:number, maxIterations?:number}} [options]
   * @returns {{matrix, residuals:number[], rmse:number, maxIndex:number, pointCount:number,
   *            inliers:boolean[], inlierCount:number, outlierCount:number,
   *            tolerancePx:number, candidatesEvaluated:number, inlierSse:number,
   *            leaveOneOut:object}}
   * @throws {FitError} code ∈ TOO_FEW_POINTS | BAD_INPUT |
   *                    NO_NONCOLLINEAR_TRIPLET | NO_NONCOLLINEAR_INLIERS
   */
  function fitAffineRobust(src, dst, options) {
    validate(src, dst);
    const n = src.length;
    const opts = options || {};
    const maxCandidates = Math.max(1, Math.floor(opts.maxCandidates || MAX_CANDIDATES));
    const maxIterations = Math.max(1, Math.floor(opts.maxIterations || MAX_REFINE_ITERS));
    const tolerancePx = opts.tolerancePx === undefined ? DEFAULT_TOLERANCE_PX : opts.tolerancePx;
    if (!Number.isFinite(tolerancePx) || tolerancePx < 0) {
      throw new FitError('BAD_INPUT', '参考图像素容差必须为不小于 0 的有限数字');
    }
    // 与坐标尺度相关的浮点容差（tol=0 时三点自身仍可判为内点）
    let scaleMax = 1;
    for (const p of dst) {
      const m = Math.max(Math.abs(p.x), Math.abs(p.y));
      if (m > scaleMax) scaleMax = m;
    }
    const eps = 1e-9 * scaleMax;
    const isInlier = (r) => r <= tolerancePx + eps;

    // ---- 有上限地枚举非共线三点候选（固定字典序，保证结果可复现）----
    let best = null; // { matrix, inliers:boolean[], count, sse }
    let evaluated = 0;
    let scanned = 0;
    let triSrcDegenerate = false;
    let triDstDegenerate = false;
    let hitScanLimit = false;

    outer:
    for (let i = 0; i < n - 2; i++) {
      for (let j = i + 1; j < n - 1; j++) {
        for (let k = j + 1; k < n; k++) {
          if (scanned >= maxCandidates) { hitScanLimit = true; break outer; }
          scanned++;
          const s3 = [src[i], src[j], src[k]];
          const d3 = [dst[i], dst[j], dst[k]];
          if (collinearity(s3).degenerate) { triSrcDegenerate = true; continue; }
          if (collinearity(d3).degenerate) { triDstDegenerate = true; continue; }

          const matrix = solveAffine(s3, d3);
          const residuals = residualsOf(matrix, src, dst);
          const inliers = new Array(n).fill(false);
          let count = 0, sse = 0;
          for (let t = 0; t < n; t++) {
            if (isInlier(residuals[t])) { inliers[t] = true; count++; sse += residuals[t] * residuals[t]; }
          }
          evaluated++;
          // 内点数量优先、内点平方误差次之；枚举按固定顺序，平局保留首个
          if (!best || count > best.count || (count === best.count && sse < best.sse)) {
            best = { matrix, inliers, count, sse };
          }
        }
      }
    }

    if (!best) {
      let message;
      if (!hitScanLimit) {
        message = triSrcDegenerate && !triDstDegenerate
          ? '扫描图中找不到三对非共线控制点，无法求解仿射变换，请调整点位分布'
          : triDstDegenerate && !triSrcDegenerate
            ? '参考图中找不到三对非共线控制点，无法求解仿射变换，请调整点位分布'
            : '找不到三对非共线控制点（点位重合或近共线），无法求解仿射变换，请调整点位分布';
      } else {
        message = '已评估的控制点三元组均近共线或重合，找不到可用的非共线三点，请调整点位分布';
      }
      throw new FitError('NO_NONCOLLINEAR_TRIPLET', message);
    }

    // ---- 对内点重新最小二乘拟合，并迭代重分类至稳定 ----
    let inliers = best.inliers;
    let matrix = best.matrix;
    let refined = false;
    for (let iter = 0; iter < maxIterations; iter++) {
      const idx = [];
      for (let t = 0; t < n; t++) if (inliers[t]) idx.push(t);
      const inSrc = idx.map((t) => src[t]);
      const inDst = idx.map((t) => dst[t]);
      if (idx.length >= MIN_POINTS &&
          !collinearity(inSrc).degenerate && !collinearity(inDst).degenerate) {
        matrix = solveAffine(inSrc, inDst);
        refined = true;
      }
      // 否则保留上一步的模型（best 至少含 3 个非共线内点）

      const residuals = residualsOf(matrix, src, dst);
      const next = new Array(n).fill(false);
      let changed = false;
      for (let t = 0; t < n; t++) {
        next[t] = isInlier(residuals[t]);
        if (next[t] !== inliers[t]) changed = true;
      }
      inliers = next;
      if (!changed) break;
    }

    const inlierIdx = [];
    for (let t = 0; t < n; t++) if (inliers[t]) inlierIdx.push(t);
    if (inlierIdx.length < MIN_POINTS) {
      // 正常不会发生（容差只增不减时内点只会更稳定），保留点对并禁用导出
      throw new FitError('NO_NONCOLLINEAR_INLIERS', '容差内的非共线内点不足 3 对，无法得到稳健模型，请放宽像素容差或调整点位');
    }
    if (collinearity(inlierIdx.map((t) => src[t])).degenerate ||
        collinearity(inlierIdx.map((t) => dst[t])).degenerate) {
      throw new FitError('NO_NONCOLLINEAR_INLIERS', '容差内的内点近共线（或重合），无法求解仿射变换，请放宽像素容差或调整点位分布');
    }

    // ---- 最终模型与统计（RMSE 仅统计内点；maxIndex 为内点中残差最大者）----
    const residuals = residualsOf(matrix, src, dst);
    let sse = 0, maxIndex = -1;
    for (const t of inlierIdx) {
      sse += residuals[t] * residuals[t];
      if (maxIndex === -1 || residuals[t] > residuals[maxIndex]) maxIndex = t;
    }
    const rmse = Math.sqrt(sse / inlierIdx.length);

    // 逐点留出验证仅以当前内点为参与集；外点 status 为 excluded（不参与）
    const leaveOneOut = leaveOneOutValidate(src, dst, inliers);

    return {
      matrix,
      residuals,
      rmse,
      maxIndex,
      pointCount: n,
      inliers,
      inlierCount: inlierIdx.length,
      outlierCount: n - inlierIdx.length,
      tolerancePx,
      candidatesEvaluated: evaluated,
      inlierSse: sse,
      refined,
      leaveOneOut
    };
  }

  /** 对单个点应用六参数矩阵 */
  function applyMatrix(m, p) {
    return { x: m.a * p.x + m.b * p.y + m.c, y: m.d * p.x + m.e * p.y + m.f };
  }

  /**
   * 逐点留出验证（Leave-One-Out Cross-Validation）。
   *
   * 对每个参与点（mask[i] 为真）：移除该点后用其余参与点重新做经典最小
   * 二乘仿射拟合，再用新模型预测留出点的参考图位置，预测位置与实际参考
   * 点的欧氏距离即该点的留出验证误差。留出拟合独立完成，不使用含留出点
   * 的原拟合残差。非参与点（mask[i] 为假，如稳健模式外点）标记为不参与。
   *
   * 参与点不足 MIN_LOO_POINTS（4）个时，所有参与点均无法验证；某次留出
   * 后的点集在任一侧近共线（或重合）时，仅该点标记为无法验证。汇总只
   * 统计可验证点，零个可验证点时汇总字段为 null（不可用）。
   *
   * @param {Array<{x:number,y:number}>} src 源点（扫描图）
   * @param {Array<{x:number,y:number}>} dst 目标点（参考图）
   * @param {boolean[]} [mask] 参与集掩码，长度与点数组一致；省略时全部参与
   * @returns {{participating:number, verifiable:number, rmse:number|null,
   *            maxError:number|null, maxIndex:number|null,
   *            points:Array<{index:number, status:'verifiable'|'unverifiable'|'excluded',
   *            error:number|null, predicted:{x,y}|null, reason?:string}>}}
   * @throws {FitError} code ∈ BAD_INPUT
   */
  function leaveOneOutValidate(src, dst, mask) {
    if (!Array.isArray(src) || !Array.isArray(dst) || src.length !== dst.length) {
      throw new FitError('BAD_INPUT', '源点与目标点必须为等长数组');
    }
    const n = src.length;
    for (const p of src.concat(dst)) {
      if (!p || !Number.isFinite(p.x) || !Number.isFinite(p.y)) {
        throw new FitError('BAD_INPUT', '控制点坐标必须为有限数字');
      }
    }
    const inMask = mask || new Array(n).fill(true);
    if (!Array.isArray(inMask) || inMask.length !== n) {
      throw new FitError('BAD_INPUT', '参与集掩码必须与控制点等长');
    }

    const partIdx = [];
    for (let i = 0; i < n; i++) if (inMask[i]) partIdx.push(i);

    const points = [];
    for (let i = 0; i < n; i++) {
      if (!inMask[i]) {
        points.push({ index: i, status: 'excluded', error: null, predicted: null, reason: 'outlier' });
        continue;
      }
      // 留出后剩余参与点不足 3 个 → 无法重拟合
      if (partIdx.length < MIN_LOO_POINTS) {
        points.push({ index: i, status: 'unverifiable', error: null, predicted: null, reason: 'too_few_points' });
        continue;
      }
      const restSrc = [];
      const restDst = [];
      for (const t of partIdx) {
        if (t !== i) { restSrc.push(src[t]); restDst.push(dst[t]); }
      }
      // 留出后的点集近共线（或重合）→ 该点无法验证，不用原残差充数
      if (collinearity(restSrc).degenerate || collinearity(restDst).degenerate) {
        points.push({ index: i, status: 'unverifiable', error: null, predicted: null, reason: 'collinear' });
        continue;
      }
      const m = solveAffine(restSrc, restDst);
      const pred = applyMatrix(m, src[i]);
      points.push({
        index: i,
        status: 'verifiable',
        error: Math.hypot(pred.x - dst[i].x, pred.y - dst[i].y),
        predicted: pred
      });
    }

    let verifiable = 0, sse = 0, maxError = null, maxIndex = null;
    for (const p of points) {
      if (p.status !== 'verifiable') continue;
      verifiable++;
      sse += p.error * p.error;
      if (maxError === null || p.error > maxError) { maxError = p.error; maxIndex = p.index; }
    }

    return {
      participating: partIdx.length,
      verifiable,
      rmse: verifiable > 0 ? Math.sqrt(sse / verifiable) : null,
      maxError,
      maxIndex,
      points
    };
  }

  return {
    fitAffine,
    fitAffineRobust,
    leaveOneOutValidate,
    applyMatrix,
    collinearity,
    FitError,
    MIN_POINTS,
    MIN_LOO_POINTS,
    COLLINEAR_RATIO,
    MAX_CANDIDATES,
    MAX_REFINE_ITERS,
    DEFAULT_TOLERANCE_PX
  };
});
