#ifndef COSMETICS                                                                                                                                                                                         
#define COSMETICS    

int color_H = kYellow-6;
int color_N = kCyan-5;
int color_P = kBlue+4;

int color_ND3 = color_N;
int color_C = color_P;
int color_NH3 = color_H;

int color_accent = 801;
int color_accent2 = 38;

void style_PbPt_plots()
{
  // Create a new TStyle object
  TStyle *PbPtStyle = new TStyle("PbPtStyle", "PbPt Style");
  
  PbPtStyle->SetOptTitle(0); // Turn off titles
  PbPtStyle->SetFrameLineWidth(0); //Remove box around the frame
  PbPtStyle->SetFrameBorderMode(0);
  PbPtStyle->SetOptStat("e"); //No stat box
  PbPtStyle->SetTitleSize(0.07, "X"); // Set title size for X axis
  PbPtStyle->SetTitleSize(0.07, "Y"); // Set title size for Y axis
  PbPtStyle->SetTitleOffset(0.8, "X"); // Set title offset for X axis
  PbPtStyle->SetTitleOffset(0.8, "Y"); // Set title offset for Y axis
  
  PbPtStyle->SetLabelSize(0.06, "X"); // Set label size for X axis
  PbPtStyle->SetLabelSize(0.06, "Y"); // Set label size for Y axis
  PbPtStyle->SetNdivisions(7, "X");
  PbPtStyle->SetNdivisions(710, "Y");
  
  PbPtStyle->SetPadLeftMargin(0.15);
  PbPtStyle->SetPadBottomMargin(0.15);
  PbPtStyle->SetPadBorderSize(0);
  PbPtStyle->SetPadBorderMode(0);
  
  PbPtStyle->SetLegendBorderSize(-1); 
  PbPtStyle->SetLegendFont(42); 
  
  // Set the default style to the PbPt style
  gROOT->SetStyle("PbPtStyle");
  
  // Apply the PbPt style to the current canvas
  gROOT->ForceStyle();
}

void style_large_plots()
{
  // Create a new TStyle object
  TStyle *largeStyle = new TStyle("largeStyle", "large Style");
  
  largeStyle->SetOptTitle(0); // Turn off titles
  largeStyle->SetFrameLineWidth(0); //Remove box around the frame
  largeStyle->SetFrameBorderMode(0);
  largeStyle->SetOptStat("e"); //No stat box
  largeStyle->SetTitleSize(0.04, "X"); // Set title size for X axis
  largeStyle->SetTitleSize(0.04, "Y"); // Set title size for Y axis
  largeStyle->SetTitleOffset(1.0, "X"); // Set title offset for X axis
  largeStyle->SetTitleOffset(1.2, "Y"); // Set title offset for Y axis
  
  largeStyle->SetLabelSize(0.04, "X"); // Set label size for X axis
  largeStyle->SetLabelSize(0.04, "Y"); // Set label size for Y axis
  largeStyle->SetNdivisions(7, "X");
  largeStyle->SetNdivisions(710, "Y");
  
  largeStyle->SetPadLeftMargin(0.15);
  largeStyle->SetPadBottomMargin(0.15);
  largeStyle->SetPadBorderSize(0);
  largeStyle->SetPadBorderMode(0);
  
  largeStyle->SetLegendBorderSize(-1); 
  largeStyle->SetLegendFont(42); 

  largeStyle->SetMarkerSize(1.5);
  largeStyle->SetMarkerStyle(21);
  
  // Set the default style to the large style
  gROOT->SetStyle("largeStyle");
  
  // Apply the large style to the current canvas
  gROOT->ForceStyle();
}


#endif
