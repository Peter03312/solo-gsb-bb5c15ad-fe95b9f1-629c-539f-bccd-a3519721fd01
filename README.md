# 扫描印刷稿套色偏移校准工具

用于校准扫描印刷稿套色（Registration）偏移的 Web 应用。上传参考图与扫描图，成对标注控制点，系统求出扫描图 → 参考图的六参数仿射变换，实时显示各点残差与**逐点留出验证误差**，叠加预览校准效果，并可导出矩阵与残差 JSON。

拟合支持两种可切换的模式：

- **经典最小二乘**：全部控制点参与，最小化残差平方和（≥ 3 对）。
- **稳健模式（默认关闭，可随时开关）**：按**参考图像素容差**区分内点/外点。算法按固定顺序枚举有上限的非共线三点候选，以**内点数量优先、内点平方误差次之**选出最佳模型，再对全部内点重新最小二乘拟合并迭代至内点集合稳定。枚举上限保证计算量有界，点对增多时拖动不卡顿。外点在双图标记、残差表与导出文件中均单独标识，**不计入最终 RMSE**。

**逐点留出验证（LOOCV）**：对每个参与拟合的控制点，移除该点后用其余参与点重新做一次最小二乘仿射拟合，再预测留出点的参考图位置，预测位置与实际参考点的欧氏距离即该点的留出验证误差。该误差衡量模型对未参与该次拟合点位的预测能力，比普通（含自身的）拟合残差更能反映真实校准精度：

- 稳健模式**仅以当前内点**为参与集；外点在验证列显示「不参与」。
- 参与点不足 **4 个**（留出后不足 3 个），或某次留出后的点集近共线/重合时，对应点标记「无法验证」，**不会用原拟合残差充数**。
- 汇总（验证 RMSE、最大验证误差、可验证点数）只统计可验证点；零个可验证点时显示「不可用」。
- 点位增删/拖动、模式切换、容差调整后，验证误差列、最差验证点提示与汇总均**实时更新**；`/api/fit` 返回与页面完全一致的 `leaveOneOut` 结果（共用 `public/js/affine.js`）。

- 后端：Node.js + Express（静态托管 + `/api/fit` 拟合接口）
- 前端：原生 Canvas（无任何前端框架/构建步骤）
- 拟合计算在浏览器本地完成，图像不会上传服务器

## 快速开始

### 方式一：Docker（单条命令，含依赖安装与启动）

```bash
docker compose up --build
```

启动后在浏览器打开：**http://localhost:3000**

### 方式二：本地 Node.js（≥ 18）

```bash
npm install && npm start
```

### 访问地址

启动后在浏览器打开：**http://localhost:3000**（端口可通过环境变量 `PORT` 修改；Docker 端口映射同样为 `3000:3000`）

## 使用说明

1. **上传图像**：分别选择参考图与扫描图（或点击「载入示例图像」快速体验，示例内置已知偏移）。
2. **标注控制点**：在任意一幅图上点击开始一对点，再在另一幅图上点击对应位置完成配对；至少 **3 对**才能拟合，建议 **4 对以上**且点位尽量分散——留出验证至少需要 4 个参与点（3 对时为精确解、拟合残差恒为 0，且无法做留出验证）。
3. **微调与删除**：拖动标记微调位置；右键标记或点击列表中的 ✕ 删除；`Esc` 取消未完成的半对点。任何点位变化都会**立即重算**拟合残差与留出验证。
4. **开启稳健模式**：勾选工具栏「**稳健模式**」，并在「**像素容差**」中填写内点允许的残差上限（参考图像素，默认 2 px，可填 0）。容差调大/调小、点位编辑或删除都会**同步重算**（内点身份变化时，留出验证的参与集随之更新），无需重新打点。
5. **查看结果**：右侧面板显示六参数矩阵、RMSE 与留出验证汇总（验证 RMSE / 最大验证误差 / 可验证点数）；残差表含**残差**与**留出验证 (px)** 两列：内点最大拟合残差橙色标出，**留出验证误差最大的可验证点青色加粗**并在状态条提示「最差验证点」；**外点**红色虚线、红底「外点」标签且验证列显示「不参与」；「无法验证」灰色标注原因（参与点不足或留出后近共线）。稳健模式的 RMSE 与验证汇总仅统计内点。控制点近共线（或重合）导致找不到至少三对非共线内点时，**点对会保留**，界面说明原因并**禁用导出**，调整点位或放宽容差后自动恢复。
6. **叠加预览**：底部画布将扫描图按拟合矩阵变换后与参考图叠加，用透明度滑杆对比校准效果。
7. **导出**：可下载六参数矩阵 JSON（含行主序 3×3 形式与留出验证汇总）与逐点残差 JSON（每点含拟合残差、内/外点身份、留出验证状态与误差、是否最差验证点）。叠加预览与两种导出使用**同一最终矩阵**，导出文件记录模式、容差、每对点身份及留出验证结果。

## 拟合模型

求解六参数仿射变换：

```
u = a·x + b·y + c
v = d·x + e·y + f        (x, y)：扫描图坐标 → (u, v)：参考图坐标
```

实现上先对坐标中心化，正规方程退化为两个 2×2 线性方程组，避免大坐标（高分辨率扫描）下的数值病态。近共线判定基于点集协方差矩阵最小/最大特征值之比（阈值 `1e-4`，见 `public/js/affine.js` 中 `COLLINEAR_RATIO`）。

稳健模式（`fitAffineRobust`）流程：

1. 按 `(i, j, k)` 字典序扫描至多 `MAX_CANDIDATES`（256）个三元组，跳过任一侧近共线/重合者；
2. 每个候选由三点求出精确仿射模型，统计残差 ≤ `tolerancePx` 的内点；
3. 以内点数量最多、次以内点平方误差和最小选择模型（平局保留固定顺序中的首个，结果可复现）；
4. 对内点重新做最小二乘拟合，并用新模型重新分类，迭代（至多 `MAX_REFINE_ITERS` 次）至内点集合不变；
5. 输出全部点的残差，但 RMSE 与 `maxResidualIndex` 只基于内点；内点不足 3 对或近共线时报错。

逐点留出验证（`leaveOneOutValidate`，两种拟合模式都会随结果返回 `leaveOneOut`）：

1. 参与集：经典模式为全部控制点；稳健模式为最终内点集合（外点 `status: "excluded"`、`reason: "outlier"`）；
2. 参与点少于 `MIN_LOO_POINTS`（4）时，参与点全部 `status: "unverifiable"`、`reason: "too_few_points"`；
3. 否则对每个参与点，用其余参与点单独做一次经典最小二乘拟合（不做 RANSAC、不使用含该点的原模型），预测该点参考图位置并记录欧氏误差；留出后点集任一侧近共线/重合时仅该点 `unverifiable`、`reason: "collinear"`；
4. 汇总仅统计 `verifiable` 点：`rmse`（验证误差均方根）、`maxError`、`maxIndex`；零个可验证点时三者为 `null`。

## API

拟合逻辑（`public/js/affine.js`）在浏览器与 Node 间共用，服务器同时暴露接口供程序化调用：

```bash
# 经典最小二乘
curl -X POST http://localhost:3000/api/fit \
  -H 'Content-Type: application/json' \
  -d '{"src":[{"x":0,"y":0},{"x":100,"y":0},{"x":0,"y":100}],
       "dst":[{"x":15,"y":-20},{"x":117,"y":-19},{"x":14,"y":78}]}'

# 稳健模式（容差 2px）
curl -X POST http://localhost:3000/api/fit \
  -H 'Content-Type: application/json' \
  -d '{"robust":true,"tolerancePx":2,
       "src":[{"x":0,"y":0},{"x":100,"y":0},{"x":0,"y":100},{"x":500,"y":500}],
       "dst":[{"x":15,"y":-20},{"x":117,"y":-19},{"x":14,"y":78},{"x":0,"y":0}]}'
```

| 接口 | 说明 |
| --- | --- |
| `GET /api/health` | 健康检查 |
| `POST /api/fit` | 仿射拟合。请求体 `{src, dst, robust?, tolerancePx?}`，`src/dst` 为等长控制点数组（≥ 3）。成功返回 `{ok, mode, matrix, residuals, rmse, maxIndex, pointCount, inliers, leaveOneOut, ...}`；`leaveOneOut` 为逐点留出验证结果（`participating/verifiable/rmse/maxError/maxIndex/points`，不可用时汇总字段为 `null`）。稳健模式额外含 `inlierCount/outlierCount/tolerancePx/candidatesEvaluated`。参数非法 `400`，点数不足/近共线/无可用非共线三点 `422`。也可直接调用 `POST` 后使用 `Affine.leaveOneOutValidate(src, dst, mask?)` 复用同一验证逻辑 |

## 导出文件格式

**affine-matrix.json**

```json
{
  "mode": "robust",
  "matrix": { "a": 1.0, "b": 0.0, "c": 0.0, "d": 0.0, "e": 1.0, "f": 0.0 },
  "matrixRowMajor": [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
  "rmse": 0.123,
  "rmseScope": "inliers",
  "leaveOneOut": {
    "participating": 4,
    "verifiable": 4,
    "rmse": 0.187,
    "maxError": 0.312,
    "maxErrorIndex": 3,
    "points": [
      { "index": 1, "status": "verifiable", "leaveOneOutErrorPx": 0.11 },
      { "index": 5, "status": "excluded", "leaveOneOutErrorPx": null, "reason": "outlier" }
    ]
  },
  "pointCount": 5,
  "tolerancePx": 2,
  "inlierCount": 4,
  "outlierCount": 1
}
```

**residuals.json**：`{ mode, rmse, rmseScope, tolerancePx?, inlierCount, outlierCount, maxResidualIndex, leaveOneOut, points: [{ id, scan, ref, residualPx, isInlier, isMaxResidual, leaveOneOutStatus, leaveOneOutErrorPx, leaveOneOutReason?, isWorstLeaveOneOut }] }`。`leaveOneOutStatus` 取值 `verifiable`（有验证误差）、`unverifiable`（`too_few_points`/`collinear`，误差为 `null`）或 `excluded`（稳健模式外点，误差为 `null`）。外点 `isInlier` 为 `false`，其拟合残差仍列出但不参与 `rmse` 与留出验证。顶层 `leaveOneOut.rmse/maxError/maxErrorIndex` 在零个可验证点时为 `null`。

## 项目结构

```
├── server.js              # Express 入口：静态托管 + /api/fit（支持 robust 参数，返回 leaveOneOut）
├── public/
│   ├── index.html         # 页面（稳健模式开关 + 像素容差输入 + 留出验证列）
│   ├── css/style.css
│   └── js/
│       ├── affine.js      # 经典/稳健仿射拟合 + 逐点留出验证（浏览器/Node 共用）
│       └── app.js         # 画布交互、实时重算、内外点标识、留出验证展示、叠加预览、导出
├── test/
│   ├── affine.test.js     # 拟合数学单元测试（含稳健模式与留出验证）
│   └── ui.smoke.js        # 前端全流程冒烟测试（DOM 桩，含稳健模式与留出验证）
├── Dockerfile
├── docker-compose.yml
└── package.json
```

## 测试

```bash
npm test
```

覆盖：无噪声精确还原、含噪最小二乘逼近、点数不足/共线/重合/非法输入拒绝、大坐标数值稳定性；稳健模式的外点剔除、内点重拟合、容差影响、固定顺序与候选上限、无非共线三点拒绝；逐点留出验证的无噪声精确预测、3 点/近共线「无法验证」且汇总为 `null`（不以原残差充数）、LOO RMSE ≥ 拟合 RMSE 的 PRESS 性质、稳健模式外点「不参与」、独立 API 的参与掩码与非法入参；以及前端完整交互流程（打点→拟合与验证列/汇总/最差点提示→导出含验证结果→删点后不可验证→共线拒绝→稳健模式切换、外点不参与验证与容差调整）。

## 说明

- 更换任一图像会清空已有控制点（坐标随图像失效）。
- 画布按图像原始分辨率渲染、CSS 自适应缩放，点击坐标按缩放比换算回图像像素，保证取点精度。
- 稳健模式只改变拟合与标识，不会删除任何点对；关闭开关即回到经典最小二乘。
- 留出验证是「诚实」的预测误差：每次重拟合都不含被预测点本身；3 对控制点（或稳健模式仅 3 个内点）时没有统计自由度，界面明确显示「无法验证 / 不可用」，而非用拟合残差代替。
