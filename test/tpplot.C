{
  TCanvas MyCanvas("MyCanvas", "New Geometry");
  MyCanvas->Divide(1,2);
  MyCanvas->cd(1);
  SHNtuple->Draw("simtrk_simhit.gy:simtrk_simhit.gx", "simtrk_simhit.subid == 1 || simtrk_simhit..subid == 3 || simtrk_simhit.subid == 5");
  htemp->SetYTitle("Y (cm)");htemp->SetXTitle("X (cm)");htemp->SetTitle("Tracker hits |z|<30 (cm)");

  MyCanvas->cd(2);
  SHNtuple->Draw("sqrt((simtrk_simhit.gy*simtrk_simhit.gy)+(simtrk_simhit.gx*simtrk_simhit.gx)):simtrk_simhit.gz");
     TH2D *htemp = (TH2D*)gPad->GetPrimitive("htemp");
     htemp->SetYTitle("R (cm)");
     htemp->SetXTitle("Z (cm)");
     htemp->SetTitle("Tracker Hits");
     TAxis *axis = htemp->GetYaxis();
     axis->SetLimits(0., 115.);
     TAxis *xaxis = htemp->GetXaxis();
     xaxis->SetLimits(-320., 320.);
     MyCanvas_2->RedrawAxis();

     TLine l0= TLine(0.0, 0.0,   0.00, 115.0); l0.SetLineColor(2);
     TLine n1= TLine(0.0, 0.0, 135.12, 115.0); n1.SetLineColor(2);
     TLine l15=TLine(0.0, 0.0, 244.83, 115.0);l15.SetLineColor(2);
     TLine n2= TLine(0.0, 0.0, 320.00, 88.20); n2.SetLineColor(2);
     TLine l25=TLine(0.0, 0.0, 320.00, 52.92);l25.SetLineColor(2);
     TLine m1= TLine(0.0, 0.0, -135.12, 115.0); m1.SetLineColor(2);
     TLine m15=TLine(0.0, 0.0, -244.83, 115.0);m15.SetLineColor(2);
     TLine m2= TLine(0.0, 0.0, -320.00, 88.20); m2.SetLineColor(2);
     TLine m25=TLine(0.0, 0.0, -320.00, 52.92);m25.SetLineColor(2);
     l0->Draw("Same");
     n1->Draw("Same"); l15->Draw("Same"); n2->Draw("Same"); l25->Draw("Same");
     m1->Draw("Same"); m15->Draw("Same"); m2->Draw("Same"); m25->Draw("Same");
     TLatex l; l.SetTextAlign(12); l.SetTextSize(0.04); l.SetTextColor(1);
     l.DrawLatex(5,110,"#eta = 0");
     l.DrawLatex(137,110,"#eta = 1");
     l.DrawLatex(247,110,"#eta = 1.5");
     l.DrawLatex(285,90,"#eta = 2");
     l.DrawLatex(285,43,"#eta = 2.5");
}

