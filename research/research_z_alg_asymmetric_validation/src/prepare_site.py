import json
from accounting import P,write

def main():
    figures=[
        ('01_geometry','非对称固定标签几何与参考材料','新六实材料、两照明、单复接收合同；不是原P1。'),
        ('02_material_response','实际材料位移与完整响应','所有节点、失败和1e-10门均保留；非零材料变化不等于decoder gauge。'),
        ('03_paired_coverage','两图册覆盖逐目标配对','固定网格与径向节点资格相同；不是全域连通分量证书。'),
        ('04_cost_windows','两个成本窗口及强普通基线','原记录窗口与统一完整trace窗口分列；总13.6%未过20%，28.6%仅子项。'),
        ('05_branch_search','独立满秩分支搜索','33启动31合格仅1根；两失败保留，往返路径不称monodromy。'),
        ('06_probability','分支权重与局部coarea积分','旧对称全域小模型与新局部核是不同对象；没有全域遗漏质量界。'),
        ('07_representation','完整换基与截断支持','完整DFT/Haar往返不创造材料；native截断违法box，其余截断也改变材料支持。'),
        ('08_certificate_scope','小区间证书的实际覆盖范围','只认证binary64系数模型的小box，不能覆盖更大的0.05数值积分片或连续PDE。')]
    report=P/'RESEARCH_REPORT_ZH.md';text=report.read_text()
    insertions={'## 4.':figures[:1],'## 6.':figures[1:3],'## 7.':figures[3:4],'## 8.':figures[4:5],'## 9.':figures[5:6],'## 10.':figures[6:]}
    for prefix,items in insertions.items():
        images='\n\n'.join('!['+title+'：'+caption+'](figures/'+stem+'.png)' for stem,title,caption in items)
        if not all('](figures/'+stem+'.png)' in text for stem,_,_ in items):text=text.replace(prefix,images+'\n\n'+prefix,1)
    report.write_text(text)
    r=json.loads((P/'RESULTS_SUMMARY.json').read_text())
    r['site']=dict(title='Z-ALG-E1｜非对称Maxwell材料图册与条件积分',subtitle='真实6-cell与4-cell · 完整实验封账 · 当前电流法总成本门未通过',
        intro_markdown='本轮确实构造了不同材料、同完整响应的非对称局部图册；直接材料法已经做到相同覆盖。电流法省了RHS，但公平总成本只省13.6%，不足预定20%。图册存在、概率正确与算法划算，分别裁决。',
        highlights=[dict(label='保守合并预算',value='21,279 / 30,000',note='15,115实际RHS＋6,164闭式等价',identity='material_geometry'),dict(label='相同局部网格覆盖',value='147 / 147',note='两方法均通过；不是147独立parent',identity='decoder_quotient'),dict(label='公平总成本门',value='13.6% < 20%',note='求解＋验收子项28.6%另列',identity='task_scope')],
        figures=[dict(path='figures/'+stem+'.png',title=title,caption=caption,sources=['figures/source_data/'+stem+'.csv']) for stem,title,caption in figures],
        required_figures=['figures/'+x[0]+'.png' for x in figures],
        appendices=[dict(path=n) for n in ['THEORY_TO_EXPERIMENT_AUDIT.md','MATHEMATICAL_PROOFS.md','RESOURCE_ACCOUNTING.md','EVIDENCE_PACKAGE.md','REPRODUCE.md','audit/prior_art/PUBLIC_CITATION_REVIEW.md','audit/history_contract/AUDIT_ZH.md','audit/solver_review/REVIEW.md','audit/results_independent/final_refresh/FINAL_REVIEW_ZH.md','audit/results_independent/final_refresh/CLAIM_REVIEW_ZH.md']],
        downloads=[dict(path='https://github.com/migodam/b1-probabilistic-maxwell-imaging/releases/tag/z-alg-asymmetric-validation-20261029',label='完整实验原始证据与GPT Pro阅读包')])
    write(P/'RESULTS_SUMMARY.json',r)
if __name__=='__main__':main()
