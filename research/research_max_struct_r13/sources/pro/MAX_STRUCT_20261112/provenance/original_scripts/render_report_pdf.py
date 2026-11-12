"""Create a self-contained Chinese PDF with ReportLab and rendered equations.

Pandoc supplies a structural AST; no browser, network, or LaTeX installation is
needed. Formula PNGs are 600 dpi, while report text and tables remain PDF text.
"""
from pathlib import Path
import json, re, html, hashlib, subprocess
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.cidfonts import UnicodeCIDFont
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.utils import ImageReader
from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph,
    Spacer, PageBreak, Table, TableStyle, Image, KeepTogether, Flowable)
from reportlab.platypus.tableofcontents import TableOfContents
import matplotlib
matplotlib.use('Agg')
from matplotlib.mathtext import math_to_image
from matplotlib.font_manager import FontProperties
from PIL import Image as PILImage

ROOT=Path(__file__).resolve().parent.parent
BASE=ROOT/'output/MAX_STRUCT_20261112'
TMP=ROOT/'tmp/pdfs';TMP.mkdir(parents=True,exist_ok=True)
MATH=TMP/'math';MATH.mkdir(exist_ok=True)
OUT=ROOT/'output/pdf/MAX_STRUCT_Research_20261112.pdf';OUT.parent.mkdir(parents=True,exist_ok=True)
ast=json.loads(subprocess.run(['pandoc',str(BASE/'report.md'),
    '--from=markdown+tex_math_dollars+tex_math_single_backslash','--to=json'],
    capture_output=True,text=True,check=True).stdout)
pdfmetrics.registerFont(UnicodeCIDFont('STSong-Light'))
pdfmetrics.registerFontFamily('STSong-Light',normal='STSong-Light',bold='STSong-Light',italic='STSong-Light',boldItalic='STSong-Light')
INK=colors.HexColor('#172b3b');BLUE=colors.HexColor('#173f53');GRAY=colors.HexColor('#586e7d')
FONT='STSong-Light';PW,PH=A4;MARGIN=49;WIDTH=PW-2*MARGIN
styles={
 'body':ParagraphStyle('body',fontName=FONT,fontSize=10.3,leading=17.2,textColor=INK,
       spaceAfter=7.5,wordWrap='CJK',allowWidows=0,allowOrphans=0),
 'h1':ParagraphStyle('h1',fontName=FONT,fontSize=17,leading=24,textColor=BLUE,
       spaceBefore=17,spaceAfter=11,keepWithNext=True,wordWrap='CJK'),
 'h2':ParagraphStyle('h2',fontName=FONT,fontSize=12.5,leading=20,textColor=BLUE,
       spaceBefore=12,spaceAfter=7,keepWithNext=True,wordWrap='CJK'),
 'h3':ParagraphStyle('h3',fontName=FONT,fontSize=11.3,leading=18,textColor=BLUE,
       spaceBefore=9,spaceAfter=6,keepWithNext=True,wordWrap='CJK'),
 'caption':ParagraphStyle('caption',fontName=FONT,fontSize=8.5,leading=13.2,textColor=GRAY,
       spaceBefore=4,spaceAfter=10,wordWrap='CJK'),
 'cell':ParagraphStyle('cell',fontName=FONT,fontSize=8.6,leading=13.1,textColor=INK,wordWrap='CJK'),
 'thead':ParagraphStyle('thead',fontName=FONT,fontSize=8.8,leading=13.4,textColor=BLUE,wordWrap='CJK'),
 'title':ParagraphStyle('title',fontName=FONT,fontSize=29,leading=40,textColor=BLUE,spaceAfter=19,wordWrap='CJK'),
 'subtitle':ParagraphStyle('subtitle',fontName=FONT,fontSize=14,leading=23,textColor=GRAY,spaceAfter=18,wordWrap='CJK'),
 'small':ParagraphStyle('small',fontName=FONT,fontSize=8.6,leading=14,textColor=GRAY,spaceAfter=8,wordWrap='CJK'),
}

def normalized_tex(s):
    s=re.sub(r'\s+',' ',s).strip()
    boxed=s.startswith(r'\boxed{') and s.endswith('}')
    if boxed:s=s[len(r'\boxed{'):-1]
    s=s.replace(r'\tfrac',r'\frac')
    s=re.sub(r'\\frac\s*([0-9])\s*([0-9])',r'\\frac{\1}{\2}',s)
    for old,new in [('le','leq'),('ge','geq'),('ne','neq')]:
        s=re.sub(r'\\'+old+r'\b',r'\\'+new,s)
    s=re.sub(r'\\(mathcal|mathbb|mathscr)\s+([A-Za-z])',r'\\\1{\2}',s)
    return s,boxed

math_records=[]
def equation_asset(tex,fontsize=13.4):
    norm,boxed=normalized_tex(tex)
    key=hashlib.sha256((str(fontsize)+norm).encode()).hexdigest()[:20]
    p=MATH/(key+'.png')
    if not p.exists():
        math_to_image('$'+norm+'$',p,prop=FontProperties(size=fontsize),dpi=600,format='png',color='#173f53')
    with PILImage.open(p) as im:w,h=im.size
    math_records.append(dict(original=tex,normalized=norm,file=p.name,width_pt=w*72/600,height_pt=h*72/600))
    return p,w*72/600,h*72/600,boxed

class Equation(Flowable):
    def __init__(self,tex):
        super().__init__();self.path,self.iw,self.ih,self.boxed=equation_asset(tex)
    def wrap(self,aw,ah):
        self.width=aw;sc=min(1.,(aw-20)/self.iw)
        self.dw=self.iw*sc;self.dh=self.ih*sc;self.height=self.dh+22
        return self.width,self.height
    def draw(self):
        if self.boxed:
            self.canv.setFillColor(colors.HexColor('#edf4f6'));self.canv.setStrokeColor(colors.HexColor('#b6cdd8'))
            self.canv.roundRect((self.width-self.dw)/2-10,4,self.dw+20,self.dh+14,4,stroke=1,fill=1)
        self.canv.drawImage(str(self.path),(self.width-self.dw)/2,11,width=self.dw,height=self.dh,mask='auto')

def esc(s):return html.escape(s.replace('–','-').replace('—','-').replace('\u2011','-')).replace('²','<super>2</super>').replace('³','<super>3</super>').replace('¹','<super>1</super>')
def inlines(items):
    out=[]
    for x in items:
        t=x['t'];c=x.get('c')
        if t=='Str':out.append(esc(c))
        elif t in ['Space','SoftBreak']:out.append(' ')
        elif t=='LineBreak':out.append('<br/>')
        elif t=='Strong':out.append('<font color="#173f53">'+inlines(c)+'</font>')
        elif t=='Emph':out.append(inlines(c))
        elif t=='Code':out.append('<font face="Courier" size="8.5">'+esc(c[1])+'</font>')
        elif t=='Link':out.append('<a href="'+html.escape(c[2][0],quote=True)+'" color="#176385">'+inlines(c[1])+'</a>')
        elif t=='Math':
            p,w,h,_=equation_asset(c[1],10.5);scale=min(1,WIDTH/max(w,1))
            out.append(f'<img src="{p}" width="{w*scale:.3f}" height="{h*scale:.3f}" valign="-2"/>')
        elif t=='Quoted':out.append('“'+inlines(c[1])+'”')
        elif t in ['Superscript','Subscript']:out.append(('<super>' if t=='Superscript' else '<sub>')+inlines(c)+('</super>' if t=='Superscript' else '</sub>'))
        elif t=='Span':out.append(inlines(c[1]))
        elif t=='RawInline':out.append(esc(c[1]))
        elif t=='Image':out.append(inlines(c[1]))
        else:raise ValueError(('Unhandled inline',t,c))
    return ''.join(out)

def plain(items):
    return re.sub('<[^>]+>','',inlines(items))

class Doc(BaseDocTemplate):
    def __init__(self,*args,**kw):
        super().__init__(*args,**kw);self.hid=0
    def beforeDocument(self):
        self.hid=0
    def afterFlowable(self,f):
        if isinstance(f,Paragraph) and f.style.name=='h1':
            key='chapter-'+str(self.hid);self.hid+=1
            self.canv.bookmarkPage(key);self.canv.addOutlineEntry(f.getPlainText(),key,0)
            self.notify('TOCEntry',(0,f.getPlainText(),self.page,key))

def header_footer(canvas,doc):
    canvas.saveState()
    if doc.page>1:
        canvas.setFont('Helvetica',7.2);canvas.setFillColor(GRAY)
        canvas.drawString(MARGIN,PH-29,'MAX-STRUCT  |  RESEARCH AUDIT & FINITE MATERIAL EVENTS')
        canvas.setStrokeColor(colors.HexColor('#c8d9e2'));canvas.setLineWidth(.45)
        canvas.line(MARGIN,PH-34,PW-MARGIN,PH-34)
    canvas.setFont('Helvetica',7.5);canvas.setFillColor(GRAY)
    canvas.drawString(MARGIN,25,'2026-11-12  |  Exposed-state mechanism evidence')
    canvas.drawRightString(PW-MARGIN,25,str(doc.page))
    canvas.restoreState()

def para(nodes,style='body',bullet=None):return Paragraph(inlines(nodes),styles[style],bulletText=bullet)
def table_cells(row,style):
    result=[]
    for cell in row[1]:
        bs=cell[4];txt=[]
        for b in bs:
            if b['t'] in ['Plain','Para']:txt.append(inlines(b['c']))
            else:raise ValueError(('Unexpected table content',b))
        result.append(Paragraph('<br/>'.join(txt),styles[style]))
    return result

def figure(image_node,caption_nodes=None):
    c=image_node['c'];p=BASE/c[2][0]
    im=Image(str(p));scale=min(WIDTH/im.imageWidth,415/im.imageHeight)
    im.drawWidth=im.imageWidth*scale;im.drawHeight=im.imageHeight*scale;im.hAlign='CENTER'
    cap=Paragraph(inlines(caption_nodes if caption_nodes is not None else c[1]),styles['caption'])
    return KeepTogether([Spacer(1,5),im,cap])

def render_blocks(blocks,depth=0):
    out=[]
    for b in blocks:
        t=b['t'];c=b.get('c')
        if t=='Header':out.append(para(c[2],'h'+str(min(3,c[0]))))
        elif t in ['Para','Plain']:
            if len(c)==1 and c[0]['t']=='Math' and c[0]['c'][0]['t']=='DisplayMath':out.append(Equation(c[0]['c'][1]))
            elif len(c)==1 and c[0]['t']=='Image':out.append(figure(c[0]))
            elif any(x['t']=='Math' for x in c):
                chunk=[]
                for x in c:
                    if x['t']=='Math':
                        if chunk:out.append(para(chunk));chunk=[]
                        out.append(Equation(x['c'][1]))
                    else:chunk.append(x)
                if chunk:out.append(para(chunk))
            else:out.append(para(c))
        elif t in ['BulletList','OrderedList']:
            items=c if t=='BulletList' else c[1]
            for i,bs in enumerate(items):
                for j,x in enumerate(bs):
                    if j==0 and x['t'] in ['Para','Plain']:
                        sty=ParagraphStyle('list',parent=styles['body'],leftIndent=13,firstLineIndent=0,bulletIndent=0,spaceAfter=5)
                        out.append(Paragraph(inlines(x['c']),sty,bulletText='•' if t=='BulletList' else f'{i+1}.'))
                    else:out.extend(render_blocks([x],depth+1))
        elif t=='Table':
            rows=[]
            for r in c[3][1]:rows.append(table_cells(r,'thead'))
            nh=len(rows)
            for body in c[4]:
                for r in body[2]+body[3]:rows.append(table_cells(r,'cell'))
            n=len(rows[0]);weights=[1]*n
            if n==2:weights=[.4,.6]
            elif n==3:weights=[.25,.375,.375]
            elif n==4:weights=[.22,.28,.25,.25]
            elif n==5:weights=[.27,.19,.18,.18,.18]
            widths=[WIDTH*w/sum(weights) for w in weights]
            tab=Table(rows,colWidths=widths,repeatRows=nh,hAlign='LEFT')
            tab.setStyle(TableStyle([('VALIGN',(0,0),(-1,-1),'TOP'),('BACKGROUND',(0,0),(-1,max(nh-1,0)),colors.HexColor('#e9f1f5')),
                ('LINEABOVE',(0,0),(-1,0),.6,colors.HexColor('#9ebac8')),
                ('LINEBELOW',(0,0),(-1,-1),.35,colors.HexColor('#c8d8df')),
                ('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),
                ('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),
                ('ROWBACKGROUNDS',(0,nh),(-1,-1),[colors.white,colors.HexColor('#f7fafb')])]))
            out.extend([Spacer(1,5),tab,Spacer(1,10)])
        elif t=='Figure':
            cap=c[1][1];caption=[]
            for x in cap:
                if x['t'] in ['Para','Plain']:caption+=x['c']
            body=c[2]
            assert len(body)==1 and body[0]['t'] in ['Plain','Para']
            out.append(figure(body[0]['c'][0],caption))
        elif t=='HorizontalRule':out.append(Spacer(1,10))
        elif t=='Div':out.extend(render_blocks(c[1],depth))
        else:raise ValueError(('Unhandled block',t,c))
    return out

doc=Doc(str(OUT),pagesize=A4,leftMargin=MARGIN,rightMargin=MARGIN,topMargin=48,bottomMargin=43,
    title='MAX-STRUCT：材料分支的可访问性与可辨识性',author='Research for migodam')
doc.addPageTemplates(PageTemplate(id='main',frames=[Frame(MARGIN,43,WIDTH,PH-48-43,id='content',leftPadding=0,rightPadding=0,topPadding=0,bottomPadding=0)],onPage=header_footer))
story=[Spacer(1,43),Paragraph('MAX-STRUCT',styles['subtitle']),Paragraph('材料分支的可访问性<br/>与可辨识性',styles['title']),
    Paragraph('R4-R12 独立研究审查、有限反射事件理论<br/>与实际 Maxwell 机制验证',styles['subtitle']),
    Spacer(1,15),Paragraph('为 migodam 编制 · 2026 年 11 月 12 日',styles['body']),Spacer(1,22)]
for title,desc in [('理论','有限有界材料事件、接收伴随场联系、近似停滞误差界；严格区分搜索与辨识。'),
                   ('证据','11 起点、44 运行、10 候选池审查；3 个真实 Maxwell 状态的冻结机制实验。'),
                   ('结论','新的事件 proposal 值得独立检验；尚未证明普遍优于 G19/DR。')]:
    story.append(Paragraph(title,styles['h2']));story.append(Paragraph(desc,styles['body']))
story.extend([Spacer(1,20),Paragraph('固定提交：ec3a37da6f3ee90beb61e753067d4f639f4118f5<br/>全部新增数值例均为已暴露对象；未宣称新盲测或全局恢复保证。',styles['small']),PageBreak()])
story.append(Paragraph('目录',styles['h1']))
toc=TableOfContents();toc.levelStyles=[ParagraphStyle('toc0',fontName=FONT,fontSize=10.5,leading=20,textColor=BLUE,spaceBefore=3,wordWrap='CJK')]
story.extend([toc,PageBreak()]);story.extend(render_blocks(ast['blocks']))
doc.multiBuild(story)
(TMP/'equation_manifest.json').write_text(json.dumps(math_records,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(dict(pdf=str(OUT),equations=len(math_records),bytes=OUT.stat().st_size),ensure_ascii=False))
