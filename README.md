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

//请注意以下的数据仅针对thyroid样本预实验。

* **阶段1：数据清洗**
    * **a. 读取GTF:** `gtf_parser.go` 负责解析 `gencode.v44.annotation.gtf.gz` 文件，计算出每个基因的准确长度。
    * **b. 加工数据:** `gct_processor.go` 利用 `gtf_parser.go` 拿到的基因长度，对 `gene_reads_v10_thyroid.gct.gz` 执行一个“两遍流式处理”：
        * **i. 第一遍:** 计算TPM标准化所需要的分母（`perSampleRPKSum`）。
        * **ii. 第二遍:** 计算 `TPM` -> `Log2(TPM+1)` 转换。
        * **iii. 两次过滤：** “低表达过滤” -> “低变异过滤”。
    * **c. 产出:** `main.go` 把这个做好的数据矩阵保存为 `clean_thyroid_matrix.csv`。
 
* **阶段2：计算相关矩阵 (Step 1: Similarity)**
    * **a. 加载数据:** `main.go` 将“阶段1”产出的 `finalMatrix`直接传递给 `RunPhase2` 函数。
    * **b. 预计算:** `phase2.go` (在 `RunPhase2` 中) 首先一次性计算所有 16,746 个基因的均值 (Mean) 和标准差 (StdDev)。
    * **c. 并行计算:** `phase2.go` (在 `correlationWorker` 中) 启动一个与CPU核心数相等的“工人池” (Worker Pool)，并行处理 1.4 亿个相关性计算任务：
        * **i. 任务分发:** 一个 `goroutine` 将基因配对任务 `[i, j]` 放入 `jobs` 管道 (channel)。
        * **ii. 并发执行:** 所有“工人” `goroutine` 从管道中抓取任务，并行计算皮尔逊相关系数 `r`。
        * **iii. 结果收集:** “工人们”将计算结果 `[i, j, r]` 放入 `results` 管道，由一个 `goroutine` 异步收集并组装成最终矩阵。
    * **d. 产出:** `main.go` 调用 `writeCorrelationMatrix` 函数，将这个 `16746 x 16746` 的相关性矩阵保存为 `correlation_matrix.csv`。
