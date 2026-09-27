"""Build the research paper PDF from saved results with embedded licensed Times New Roman."""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Image,
                                Table, TableStyle, PageBreak)
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Keep TTFs out of the repository. Extract the user's licensed times32.exe
# locally, then set VCELL_TNR_DIR to that directory before building.
FONT_DIR = os.environ.get('VCELL_TNR_DIR')
if not FONT_DIR:
    raise RuntimeError('Set VCELL_TNR_DIR to locally extracted Times New Roman TTFs')
for face, filename in [('Times-Roman','Times.TTF'), ('Times-Bold','Timesbd.TTF'),
                       ('Times-Italic','Timesi.TTF'), ('Times-BoldItalic','Timesbi.TTF')]:
    pdfmetrics.registerFont(TTFont(face, os.path.join(FONT_DIR, filename)))
pdfmetrics.registerFontFamily('Times-Roman', normal='Times-Roman', bold='Times-Bold',
                              italic='Times-Italic', boldItalic='Times-BoldItalic')

TITLE = ParagraphStyle('TitleTNR', fontName='Times-Bold', fontSize=20, leading=24, alignment=TA_CENTER, spaceAfter=18)
H1 = ParagraphStyle('H1TNR', fontName='Times-Bold', fontSize=14, leading=17, spaceBefore=14, spaceAfter=6)
H2 = ParagraphStyle('H2TNR', fontName='Times-Bold', fontSize=12, leading=15, spaceBefore=10, spaceAfter=4)
BODY = ParagraphStyle('BodyTNR', fontName='Times-Roman', fontSize=11.5, leading=15.6, alignment=TA_JUSTIFY, spaceAfter=8)
EQ = ParagraphStyle('EqTNR', fontName='Times-Italic', fontSize=11, leading=15, alignment=TA_CENTER, spaceBefore=4, spaceAfter=8)
CAP = ParagraphStyle('CapTNR', fontName='Times-Italic', fontSize=9.5, leading=12, alignment=TA_CENTER, spaceAfter=10)
REF = ParagraphStyle('RefTNR', fontName='Times-Roman', fontSize=10, leading=13, spaceAfter=4)

def P(t, s=BODY): return Paragraph(t, s)
def fig(path, w=5.6*inch, caption=''):
    from PIL import Image as PILImage
    try:
        iw, ih = PILImage.open(path).size
    except Exception:
        iw, ih = 1200, 640
    h = w * ih / iw
    return [Image(path, width=w, height=h), P(caption, CAP)]
def tbl(data, caption='', widths=None):
    t = Table(data, colWidths=widths, repeatRows=1)
    t.setStyle(TableStyle([
        ('FONTNAME', (0,0), (-1,-1), 'Times-Roman'), ('FONTSIZE', (0,0), (-1,-1), 8.5),
        ('FONTNAME', (0,0), (-1,0), 'Times-Bold'),
        ('GRID', (0,0), (-1,-1), 0.4, colors.grey),
        ('BACKGROUND', (0,0), (-1,0), colors.Color(0.9,0.9,0.95)),
        ('ROWBACKGROUNDS', (0,1), (-1,-1), [colors.white, colors.Color(0.96,0.96,0.98)]),
        ('TOPPADDING', (0,0), (-1,-1), 2), ('BOTTOMPADDING', (0,0), (-1,-1), 2)]))
    out = [t]
    if caption: out.append(P(caption, CAP))
    return out
def R(name):
    return json.load(open(os.path.join(ROOT, 'results', name)))


def main():
    import content_a, content_b, content_c, content_math, content_ext, content_v2, content_audit, content_bern, content_gate, content_calibration
    story = []
    story = content_a.story_a(story, R)
    story = content_c.story_c(story, R)
    story = content_b.story_b(story, R)
    # Insert the later completed main-body analysis before References, not after appendices.
    # A 50+ page text-body requirement cannot be met by relocating appendices.
    ref_idx = next(i for i, item in enumerate(story) if isinstance(item, Paragraph) and item.text == 'References')
    assert isinstance(story[ref_idx - 1], PageBreak)
    story[ref_idx - 1:ref_idx - 1] = content_calibration.story_calibration([], R)
    story = content_math.story_math(story, R)
    story = content_ext.story_ext(story, R)
    story = content_audit.story_audit(story, R)
    story = content_v2.story_v2(story, R)
    story = content_bern.story_bern(story, R)
    story = content_gate.story_gate(story, R)
    doc = SimpleDocTemplate(os.path.join(ROOT, 'paper', 'VC2_virtual_cell_paper.pdf'),
                            pagesize=letter,
                            leftMargin=0.9*inch, rightMargin=0.9*inch,
                            topMargin=0.9*inch, bottomMargin=0.9*inch,
                            title='VC-2: Biomass-objective sensitivity and score-lineage audit')
    from reportlab.pdfgen.canvas import Canvas
    class TNRCanvas(Canvas):
        def __init__(self, *a, **k):
            k.setdefault('initialFontName', 'Times-Roman'); super().__init__(*a, **k)
    doc.build(story, canvasmaker=TNRCanvas)
    from pypdf import PdfReader
    n = len(PdfReader(os.path.join(ROOT, 'paper', 'VC2_virtual_cell_paper.pdf')).pages)
    print('PAGES:', n)


if __name__ == '__main__':
    main()
