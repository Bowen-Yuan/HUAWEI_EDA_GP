"""Generate a readable PDF preview when a local XeLaTeX engine is unavailable."""
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import cm
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, Image, PageBreak)


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
OUT = ROOT / "adaptive_pareto_end_to_end_report_preview.pdf"


def image(path: Path, width: float) -> Image:
    item = Image(str(path)); scale = width / item.drawWidth; item.drawWidth = width; item.drawHeight *= scale
    return item


def main() -> None:
    doc = SimpleDocTemplate(str(OUT), pagesize=A4, rightMargin=1.7 * cm,
                            leftMargin=1.7 * cm, topMargin=1.6 * cm, bottomMargin=1.6 * cm)
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name="Small", parent=styles["BodyText"], fontSize=8.5, leading=11))
    styles.add(ParagraphStyle(name="H", parent=styles["Heading2"], textColor=colors.HexColor("#264653")))
    story = [Paragraph("Adaptive Pareto Soft Placement: End-to-End Experiment Report", styles["Title"]),
             Paragraph("adaptec1 numerical study - 2026-07-21", styles["Normal"]), Spacer(1, .35 * cm),
             Paragraph("PDF preview generated locally because no XeLaTeX/PDFLaTeX/Tectonic executable was found. The canonical TeX source is main.tex.", styles["Small"]),
             Spacer(1, .3 * cm), Paragraph("Executive summary", styles["H"]),
             Paragraph("The requested joint target (HPWL 55-60M and strict overflow 20-30%) was not reached. The complete initialization is balanced; the normalized HPWL rescue reaches the wirelength target but raises strict overflow; the raw-gradient soft penalty lowers HPWL further but cannot recover density from the resulting 50% state.", styles["BodyText"]), Spacer(1, .25 * cm)]
    data = [["Stage", "HPWL (M)", "Legacy O (%)", "Strict O (%)"],
            ["Initialization output", "127.623", "21.745", "30.447"],
            ["HPWL transition", "57.490", "40.694", "50.539"],
            ["Joint final (lambda max 512)", "53.998", "38.548", "49.367"]]
    table = Table(data, colWidths=[6.6 * cm, 2.6 * cm, 3.0 * cm, 2.8 * cm])
    table.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#264653")),
                               ("TEXTCOLOR", (0, 0), (-1, 0), colors.white), ("GRID", (0, 0), (-1, -1), .3, colors.HexColor("#AAAAAA")),
                               ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F3F6F6")]),
                               ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"), ("FONTSIZE", (0, 0), (-1, -1), 8),
                               ("ALIGN", (1, 1), (-1, -1), "RIGHT"), ("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("TOPPADDING", (0, 0), (-1, -1), 5), ("BOTTOMPADDING", (0, 0), (-1, -1), 5)]))
    story += [table, Spacer(1, .35 * cm), image(ROOT / "figures" / "experiment_summary.png", 17 * cm), PageBreak()]
    story += [Paragraph("Workflow and initialization", styles["H"]),
              Paragraph("The new experiment ran multilevel coarsening, a portfolio of diverse initial candidates, proxy filtering, and fine-level adaptive Pareto optimization. Its minimum-HPWL solution was subsequently rescued using the verified normalized projected HPWL subgradient alpha_k = 128 row-heights / sqrt(k), with density, lambda, trust acceptance, and short-axis intervention disabled. The final stage ramps lambda quadratically from zero to 512 for 2000 joint steps and contains no overflow-only rescue.", styles["BodyText"]), Spacer(1, .2 * cm),
              image(PROJECT / "results_e2e_init" / "adaptec1" / "20260721T034317Z" / "figures" / "05_candidate_funnel.png", 16 * cm),
              Spacer(1, .2 * cm), Paragraph("The initialization produces 127.623M HPWL and 30.447% strict overflow: a useful balanced seed, but just outside the requested strict-overflow band.", styles["Small"]), PageBreak()]
    story += [Paragraph("Convergence and evaluations", styles["H"]),
              image(PROJECT / "results_e2e_hpwl57_joint512" / "adaptec1" / "20260721T040215Z" / "figures" / "01_two_phase_convergence.png", 16 * cm),
              Spacer(1, .2 * cm), Paragraph("Short-axis ablation: standalone dominant-axis exact-oracle descent accepted only seven updates in 3000 attempts, reducing strict overflow from 51.217% to 49.325% while raising HPWL from 54.998M to 77.618M. It is locally valid but reaches a piecewise-linear plateau and cannot serve as the complete density solver.", styles["BodyText"]),
              Spacer(1, .18 * cm), Paragraph("Soft-penalty ablation: lambda=256 improves the direct 55M run only to 50.235% strict overflow. The new lambda=512 run moves 50.539% to 49.367% while HPWL falls from 57.490M to 53.998M. Therefore raw-gradient soft penalty is a gentle trade-off controller, not a density-recovery mechanism.", styles["BodyText"]), Spacer(1, .18 * cm),
              Paragraph("Earlier alternating evacuation reaches 28.939% strict overflow, but only at 693.380M HPWL; lambda continuation recovers this to 594.884M / 29.154%. Capacity-only destination selection damages net locality. The next implementation should enforce an exact HPWL budget while ranking short-side local transport by overflow decrease per estimated HPWL cost.", styles["BodyText"])]
    doc.build(story)


if __name__ == "__main__":
    main()
