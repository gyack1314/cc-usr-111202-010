# 研究中心异常计数的条件最大扫描

这是研究资料维护用的独立Python项目。`data/`保存合成记录，`research_data.io`负责严格JSON读取与SHA-256来源摘要，`research_data.domain`提供当前原始资料统计，`research_data.cli`将摘要输出为JSON。字段、时间、计数和量纲见[资料约定](docs/data-contract.md)。所有研究、中心及受试者均为虚构。

Python 3.11及以上版本可直接使用，当前没有第三方依赖。在项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q research_data tests
python -m research_data.cli inspect data/normal.json
```

`inspect`只显示已有摘要，不修改来源。文件名以`invalid-`起始的资料包含故意不合法的输入，其他文件是同结构的业务边界。数据摘要与来源指纹应同时保留，原始JSON不覆盖。

## 条件最大扫描

`research_data.scan`把“搜索过多少候选区域”纳入同一个条件null下的校正，不把事后挑中的区域当成事先指定。候选区域是各中心在固定`radii`上的图距离闭球（去重，排除空集与全体中心），候选族只由图与`radii`决定并保存SHA-256指纹，不能看`observed`后裁剪。条件null固定总异常数N，计数服从概率正比于`expected`的多项分布；每次null抽样在全中心重分配N并扫描完整候选族取最大统计量，主检验与区域校正p值都对照同一份null最大值。

```sh
# auto：状态数在 exact_state_budget 内走精确枚举，否则回退蒙特卡洛
python -m research_data.cli scan data/normal.json --outdir out/normal
# 强制蒙特卡洛（基准纯Python引擎，零第三方依赖）
python -m research_data.cli scan data/normal.json --outdir out/mc --mode monte-carlo
# 可选批量矩阵加速：与基准共用候选族与概率向量，需要 pip install -e .[numpy]
python -m research_data.cli scan data/normal.json --outdir out/np --mode monte-carlo --engine numpy
```

输出目录保存`scan_result.json`（配置摘要、种子、B、候选族与指纹、观察区域排名、主p与各区域校正p、蒙特卡洛精度说明）；蒙特卡洛模式另存全部null最大统计量`null_maxima.csv`，精确模式另存按真实概率加权的`exact_null_distribution.json`。精确模式按各计数向量的真实多项概率加权而非等权；蒙特卡洛p值为`(1+大于等于观察值次数)/(B+1)`，附带标准误差与最小分辨率`1/(B+1)`，是带抽样误差的估计，不叫精确检验。N=0时所有统计量为0、p值为1；单中心图等无有效候选区域的情形以`no_candidate_regions`单独报告。口径细节见[资料约定](docs/data-contract.md)。

可使用固定镜像容器：`docker build -t research-material .`，随后`docker run --rm research-material`查看资料；`docker run --rm research-material python -m unittest discover -s tests -v`运行现有检查。全部资料在项目内，运行不需要外部账号或服务。
