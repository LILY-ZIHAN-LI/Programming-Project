# Programming-Project
Programming Project

我打算先用一个小样本（比如说某个组织，甲状腺...）跑通这个pipeline，
然后通过将所有的组织样本跑通，合并在一起进行后续分析。

RNA seq数据：
https://www.gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression
目前先采用了Gene read counts 中thyroid的部分

annotation（用于TPM正则化）：
https://www.gencodegenes.org/human/release_44.html
Comprehensive gene annotation	CHR（基因长度注释）

*请注意以下的数据仅针对thyroid样本预实验。

* **阶段1：数据清洗**
    * **a. 读取GTF:** `gtf_parser.go` 负责解析 `gencode.v44.annotation.gtf.gz` 文件，计算出每个基因的准确长度。
    * **b. 加工数据:** `gct_processor.go` 利用 `gtf_parser.go` 拿到的基因长度，对 `gene_reads_v10_thyroid.gct.gz` 执行一个“两遍流式处理”：
        * **i. 第一遍:** 计算TPM标准化所需要的分母（`perSampleRPKSum`）。
        * **ii. 第二遍:** 计算 `TPM` -> `Log2(TPM+1)` 转换。
        * **iii. 两次过滤：** “低表达过滤” -> “低变异过滤”。
    * **c. 产出:** `main.go` 把这个做好的数据矩阵保存为 `clean_thyroid_matrix.csv`。
 
* **阶段2：计算相关矩阵**
    * **a. 加载数据:** `main.go` 将“阶段1”产出的 `finalMatrix`直接传递给 `RunPhase2` 函数。
    * **b. 预计算:** `phase2.go` (在 `RunPhase2` 中) 首先一次性计算所有 16,746 个基因的均值 (Mean) 和标准差 (StdDev)。
    * **c. 并行计算:** `phase2.go` (在 `correlationWorker` 中) 启动一个与CPU核心数相等的“工人池” (Worker Pool)，并行处理 1.4 亿个相关性计算任务：
        * **i. 任务分发:** 一个 `goroutine` 将基因配对任务 `[i, j]` 放入 `jobs` 管道 (channel)。
        * **ii. 并发执行:** 所有“工人” `goroutine` 从管道中抓取任务，并行计算皮尔逊相关系数 `r`。
        * **iii. 结果收集:** “工人们”将计算结果 `[i, j, r]` 放入 `results` 管道，由一个 `goroutine` 异步收集并组装成最终矩阵。
    * **d. 产出:** `main.go` 调用 `writeCorrelationMatrix` 函数，将这个 `16746 x 16746` 的相关性矩阵保存为 `correlation_matrix.csv`。
 
* **阶段3：邻接矩阵（Adjacency Matrix）**
    * **a. 加载数据:** `main.go` 将"阶段2"产出的 `correlationMatrix` 传递给 `CalculateAdjacencyMatrix` 函数。
    * **b. 软阈值处理:** 使用 `softPowerBeta = 6.0` 作为幂次参数，对相关性矩阵进行软阈值转换：
        * **i. 幂次转换:** 对相关性系数取绝对值后进行幂次运算：`adjacency[i][j] = |correlation[i][j]|^beta`。
        * **ii. 目的:** 强化强相关性，弱化弱相关性，使网络更符合无尺度网络特性。
    * **c. 产出:** `main.go` 调用 `writeCorrelationMatrix` 函数，将这个 `16746 x 16746` 的邻接矩阵保存为 `adjacency_matrix.csv`。

* **阶段4：拓扑重叠矩阵（TOM）**
    * **a. 加载数据:** `main.go` 将"阶段3"产出的 `adjacencyMatrix` 传递给 `CalculateTOM` 函数。
    * **b. 拓扑重叠计算:** 基于邻接矩阵计算拓扑重叠度，衡量基因对之间的共享邻居程度：
        * **i. 共享连接:** 计算基因 `i` 和基因 `j` 与其他所有基因的共享连接强度总和。
        * **ii. TOM公式:** `TOM[i][j] = (共享连接 + adjacency[i][j]) / (min(连接度_i, 连接度_j) + 1 - adjacency[i][j])`。
        * **iii. 目的:** 考虑网络拓扑结构，使相关性度量更加稳健。
    * **c. 产出:** `main.go` 调用 `writeCorrelationMatrix` 函数，将这个 `16746 x 16746` 的TOM矩阵保存为 `tom_matrix.csv`。

* **阶段5：异质性矩阵（Dissimilarity Matrix）**
    * **a. 加载数据:** `main.go` 将"阶段4"产出的 `tomMatrix` 传递给 `CalculateDissimilarity` 函数。
    * **b. 距离转换:** 将拓扑重叠度转换为距离度量：
        * **i. 转换公式:** `dissimilarity[i][j] = 1 - TOM[i][j]`。
        * **ii. 目的:** 将相似性度量转换为距离度量，为后续的层次聚类分析做准备。
    * **c. 产出:** `main.go` 调用 `writeCorrelationMatrix` 函数，将这个 `16746 x 16746` 的异质性矩阵保存为 `dissimilarity_matrix.csv`，作为最终输出供 RShiny 可视化工具进行基因模块聚类分析。


thyroid matrix:
https://drive.google.com/file/d/15FEyBlubkOzGZNARPno2PqadztDV3ztq/view?usp=drive_link
correlation matrix (Based on Pearsons):
