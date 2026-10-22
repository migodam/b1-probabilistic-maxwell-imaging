# Research P3-L PDF 构建

本目录包含独立的排版文件；脚本读取 Markdown 快照，不会修改正文。生成的 TeX、完整日志与快照写入脚本所在目录的 `latex/`，构建摘要写入 `build_summary.json`。正文中的相对图片路径以正文所在目录解析。

## 从研究包重建

在研究包目录执行（将字体目录替换为实际位置）：

```bash
python typesetting/build_report.py --source Research_P3L_Report_ZH.md --output Research_P3L_Report_ZH.pdf --font-dir /path/to/fonts
```

也支持参数名 `--input`，与 `--source` 等价。相对的输入、输出、字体路径均相对命令执行目录解释。`--passes 3` 为默认值，至少两次编译以填充目录；无需 `shell-escape`。

## 依赖与字体

需要 Python 3、Pandoc 3 与 XeLaTeX，以及常见 TeX 包 `fontspec`、`unicode-math`、`geometry`、`graphicx`、`adjustbox`、`float`、`longtable`、`booktabs`、`array`、`calc`、`etoolbox`、`needspace`、`titlesec`、`fancyhdr`、`fvextra`、`xurl` 和 `bookmark`。数学字体使用 TeX 自带的 Latin Modern Math；伪代码使用 Latin Modern Mono。**不依赖 xeCJK，不会自动安装依赖或访问网络。**

字体目录中应存在以下原文件名：

- `NotoSerifCJKsc-Regular.otf`
- `NotoSerifCJKsc-Bold.otf`

官方字体目录：<https://github.com/notofonts/noto-cjk/tree/main/Serif/OTF/SimplifiedChinese>

对应官方文件：

- <https://raw.githubusercontent.com/notofonts/noto-cjk/main/Serif/OTF/SimplifiedChinese/NotoSerifCJKsc-Regular.otf>
- <https://raw.githubusercontent.com/notofonts/noto-cjk/main/Serif/OTF/SimplifiedChinese/NotoSerifCJKsc-Bold.otf>

字体不随研究包重复分发；其授权与来源以官方仓库为准。可将下载的字体放入研究包 `fonts/`，或使用 `--font-dir` 指向已有目录。

## 排版约定与核验

A4 单栏，11 pt 中文正文，紧凑两级目录、页码和简洁页眉。技术表保留可跨页的 `longtable` 与重复表头；图按正文次序就地分页。显示公式先按自然宽度排版，仅在超出当前行宽时缩小；构建摘要记录实际缩放比例。文档标题、正文与数学内容不由构建脚本重写。

最后应检查 `build_summary.json` 中的 `warnings`、`scaled_displays` 与 `source_unchanged_during_build`，并渲染 PDF 检查实际页面。`source_unchanged_during_build=false` 表示编译期间正文发生更新，应再运行一次。`Package amsmath Warning: Foreign command \\atopwithdelims` 源于当前正文中的传统二项式写法；已检查其实际数学渲染正常。

当前已验证工具版本：Pandoc 3.1.3 与 TeX Live 2022 / XeTeX。不同 TeX 版本可有少量分页差异，正文 Markdown 和研究结果不受影响。
