// QC5graph.C
// Plot effective gas gain and trigger rate versus divider current for QC5 sectors.
//
// Run with: root -l -q 'QC5graph.C()'

#include <TCanvas.h>
#include <TAxis.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TPad.h>
#include <TSpline.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kElectronCharge = 1.6e-19; // C
constexpr double kPrimaryElectrons = 346.0;
constexpr int kPlateauWindow = 3;

constexpr Color_t kGainColor = kRed + 1;
constexpr Color_t kRateColor = kAzure + 2;
constexpr Color_t kFitColor = kGreen + 2;

struct Row {
  double imon = 0.0;
  double xoffWdt = 0.0;
  double xonWdt = 0.0;
  double ir = 0.0;
  double rate = 0.0;
  double rateError = 0.0;
  double gain = 0.0;
};

bool ParseDouble(const std::string& text, double& value) {
  try {
    std::size_t parsed = 0;
    value = std::stod(text, &parsed);
    return parsed == text.size();
  } catch (const std::exception&) {
    return false;
  }
}

// Required zero-based CSV columns: 2=Imon, 4=XOFF_WDT, 6=XON_WDT, 9=IR.
bool ReadCSV(const char* filename, std::vector<Row>& rows) {
  std::ifstream input(filename);
  if (!input) {
    std::cerr << "ERROR: Cannot open " << filename << '\n';
    return false;
  }

  std::string line;
  std::getline(input, line); // Header row.

  unsigned int lineNumber = 1;
  while (std::getline(input, line)) {
    ++lineNumber;
    if (line.empty()) continue;

    std::vector<std::string> columns;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) columns.push_back(field);

    if (columns.size() <= 9) {
      std::cerr << "WARNING: Skipping line " << lineNumber
                << " in " << filename << " (expected at least 10 columns).\n";
      continue;
    }

    Row row;
    if (!ParseDouble(columns[2], row.imon) ||
        !ParseDouble(columns[4], row.xoffWdt) ||
        !ParseDouble(columns[6], row.xonWdt) ||
        !ParseDouble(columns[9], row.ir)) {
      std::cerr << "WARNING: Skipping non-numeric line " << lineNumber
                << " in " << filename << '\n';
      continue;
    }

    row.rate = (row.xonWdt - row.xoffWdt) / 60.0;
    row.rateError = std::sqrt(std::abs(row.xonWdt + row.xoffWdt)) / 60.0;
    if (row.rate > 0.0) rows.push_back(row);
  }
  return true;
}

bool FindPlateau(const std::vector<Row>& rows, int& start, int& end,
                 double& plateauRate) {
  const int numberOfRows = static_cast<int>(rows.size());
  if (numberOfRows == 0) return false;

  const int window = std::min(kPlateauWindow, numberOfRows);
  start = 0;
  double smallestVariation = std::numeric_limits<double>::infinity();

  for (int candidate = 0; candidate <= numberOfRows - window; ++candidate) {
    double variation = 0.0;
    for (int index = candidate; index < candidate + window - 1; ++index) {
      variation += std::abs(rows[index + 1].rate - rows[index].rate);
    }
    if (variation < smallestVariation) {
      smallestVariation = variation;
      start = candidate;
    }
  }

  end = start + window - 1;
  plateauRate = 0.0;
  for (int index = start; index <= end; ++index) plateauRate += rows[index].rate;
  plateauRate /= window;
  return plateauRate > 0.0;
}

void DrawSector(const char* csvFile, const char* label) {
  std::vector<Row> rows;
  if (!ReadCSV(csvFile, rows) || rows.empty()) {
    std::cerr << "ERROR: No valid measurements in " << csvFile << '\n';
    return;
  }

  int plateauStart = 0;
  int plateauEnd = 0;
  double plateauRate = 0.0;
  if (!FindPlateau(rows, plateauStart, plateauEnd, plateauRate)) {
    std::cerr << "ERROR: Could not determine a positive rate plateau for "
              << label << '\n';
    return;
  }

  std::cout << label << ": plateau at data rows " << plateauStart << "-"
            << plateauEnd << ", rate = " << plateauRate << " Hz\n";

  const int numberOfRows = static_cast<int>(rows.size());
  std::vector<double> imon(numberOfRows), gain(numberOfRows), rate(numberOfRows);
  std::vector<double> rateError(numberOfRows), zero(numberOfRows, 0.0);
  for (int index = 0; index < numberOfRows; ++index) {
    rows[index].gain = std::abs(rows[index].ir) /
                       (kElectronCharge * kPrimaryElectrons * plateauRate);
    imon[index] = rows[index].imon;
    gain[index] = rows[index].gain;
    rate[index] = rows[index].rate;
    rateError[index] = rows[index].rateError;
  }

  const double xMin = *std::min_element(imon.begin(), imon.end());
  const double xMax = *std::max_element(imon.begin(), imon.end());
  const double xSpan = (xMax > xMin) ? xMax - xMin : 1.0;
  const double xPadding = 0.04 * xSpan;
  const double gainMin = *std::min_element(gain.begin(), gain.end());
  const double gainMax = *std::max_element(gain.begin(), gain.end());
  const double rateMin = *std::min_element(rate.begin(), rate.end());
  const double rateMax = *std::max_element(rate.begin(), rate.end());

  TCanvas* canvas = new TCanvas(Form("cQC5_%s", label), label, 1400, 800);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.14);
  canvas->SetBottomMargin(0.14);
  canvas->SetTopMargin(0.10);

  TPad* gainPad = new TPad(Form("gainPad_%s", label), "", 0, 0, 1, 1);
  gainPad->SetLeftMargin(canvas->GetLeftMargin());
  gainPad->SetRightMargin(canvas->GetRightMargin());
  gainPad->SetBottomMargin(canvas->GetBottomMargin());
  gainPad->SetTopMargin(canvas->GetTopMargin());
  gainPad->SetLogy();
  gainPad->SetTicks(1, 1);
  gainPad->SetGrid();
  gainPad->Draw();
  gainPad->cd();

  TGraphErrors* gainGraph = new TGraphErrors(
      numberOfRows, imon.data(), gain.data(), zero.data(), zero.data());
  gainGraph->SetTitle(Form("%s;Divider Current [#muA];Effective Gain", label));
  gainGraph->SetMarkerStyle(22);
  gainGraph->SetMarkerSize(1.2);
  gainGraph->SetMarkerColor(kGainColor);
  gainGraph->SetLineColor(kGainColor);
  gainGraph->GetXaxis()->SetLimits(xMin - xPadding, xMax + xPadding);
  gainGraph->GetYaxis()->SetRangeUser(std::max(gainMin * 0.7, 1e-12), gainMax * 1.3);
  gainGraph->GetYaxis()->SetTitleOffset(1.4);
  gainGraph->GetYaxis()->SetTitleColor(kGainColor);
  gainGraph->GetYaxis()->SetLabelColor(kGainColor);
  gainGraph->Draw("AP");

  TSpline3* gainSpline = new TSpline3(Form("gainSpline_%s", label), gainGraph);
  gainSpline->SetLineColor(kGainColor);
  gainSpline->SetLineWidth(2);
  gainSpline->Draw("same");

  TF1* gainFit = new TF1(Form("gainFit_%s", label), "[0]*exp([1]*x)", xMin, xMax);
  gainFit->SetParameters(std::max(gainMin, 1.0), 0.05);
  gainFit->SetLineColor(kFitColor);
  gainFit->SetLineWidth(2);
  gainGraph->Fit(gainFit, "RQ");
  gainFit->Draw("same");
  gainPad->Update();

  canvas->cd();
  TPad* ratePad = new TPad(Form("ratePad_%s", label), "", 0, 0, 1, 1);
  ratePad->SetFillStyle(4000);
  ratePad->SetFrameFillStyle(0);
  ratePad->SetLeftMargin(canvas->GetLeftMargin());
  ratePad->SetRightMargin(canvas->GetRightMargin());
  ratePad->SetBottomMargin(canvas->GetBottomMargin());
  ratePad->SetTopMargin(canvas->GetTopMargin());
  ratePad->Draw();
  ratePad->cd();

  TGraphErrors* rateGraph = new TGraphErrors(
      numberOfRows, imon.data(), rate.data(), zero.data(), rateError.data());
  rateGraph->SetMarkerStyle(20);
  rateGraph->SetMarkerSize(1.0);
  rateGraph->SetMarkerColor(kRateColor);
  rateGraph->SetLineColor(kRateColor);
  rateGraph->SetTitle("");
  rateGraph->GetXaxis()->SetLimits(xMin - xPadding, xMax + xPadding);
  rateGraph->GetYaxis()->SetRangeUser(rateMin * 0.8, rateMax * 1.2);
  rateGraph->GetYaxis()->SetLabelSize(0);
  rateGraph->GetYaxis()->SetTickSize(0);
  rateGraph->GetXaxis()->SetLabelSize(0);
  rateGraph->GetXaxis()->SetTickSize(0);
  rateGraph->Draw("AP");

  TSpline3* rateSpline = new TSpline3(Form("rateSpline_%s", label), rateGraph);
  rateSpline->SetLineColor(kRateColor);
  rateSpline->SetLineWidth(2);
  rateSpline->Draw("same");
  ratePad->Update();

  const double rightEdge = ratePad->GetUxmax();
  TGaxis* rateAxis = new TGaxis(rightEdge, rateMin * 0.8, rightEdge, rateMax * 1.2,
                                rateMin * 0.8, rateMax * 1.2, 505, "+L");
  rateAxis->SetLineColor(kRateColor);
  rateAxis->SetLabelColor(kRateColor);
  rateAxis->SetTitleColor(kRateColor);
  rateAxis->SetTitle("Rate [Hz]");
  rateAxis->SetTitleOffset(1.35);
  rateAxis->Draw();

  canvas->cd();
  TLegend* legend = new TLegend(0.55, 0.18, 0.80, 0.38);
  legend->SetBorderSize(1);
  legend->SetFillColor(kWhite);
  legend->SetTextFont(42);
  legend->SetTextSize(0.032);
  legend->AddEntry(gainGraph, "Effective gain data", "p");
  legend->AddEntry(gainFit, "Gain fit: G = A e^{#alpha I}", "l");
  legend->AddEntry(rateGraph, "Rate measurement", "p");
  legend->Draw();

  TLatex text;
  text.SetNDC();
  text.SetTextFont(42);
  text.SetTextSize(0.027);
  text.DrawLatex(0.18, 0.85, "ME0 QC5 Test");
  text.DrawLatex(0.18, 0.822, Form("Sector: %s", label));
  text.DrawLatex(0.18, 0.794, "Gas: Ar/CO_{2} (70/30)");
  text.SetTextSize(0.026);
  text.DrawLatex(0.18, 0.760, Form("Gain fit: G = %.3e #times e^{%.4f I}",
                                   gainFit->GetParameter(0), gainFit->GetParameter(1)));
  text.DrawLatex(0.18, 0.732, Form("Plateau rate: %.1f Hz  (I_{mon} %.0f-%.0f #muA)",
                                   plateauRate, rows[plateauStart].imon,
                                   rows[plateauEnd].imon));

  canvas->SaveAs(Form("%s.png", label));
  canvas->SaveAs(Form("%s.pdf", label));
  std::cout << "Saved: " << label << ".png and " << label << ".pdf\n";
}

} // namespace

void QC5graph() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridColor(kGray);

  std::cout << "QC5 gain analysis\n"
            << "Constants: e = " << kElectronCharge << " C, N = "
            << kPrimaryElectrons << ", plateau window = " << kPlateauWindow << "\n\n";
  DrawSector("sector3_2.csv", "Sector 3");
  DrawSector("sector4_2.csv", "Sector 4");
  DrawSector("sector5_2.csv", "Sector 5");
}
