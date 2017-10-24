#include <TH2F.h>
#include <TH1F.h>
#include <TFile.h>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <TCanvas.h>
#include "TStyle.h"
#include <TLine.h>
#include <TLegend.h>

// B dataset: 193833-196531
// C dataset: 198022-203742
// D dataset: 203777-208686

void lumiplotting(TFile* f)
{

	gStyle->SetOptStat(000000000);

	int bdown = 193833, bup = 196531, cdown = 198022, cup = 203742, ddown =203777, dup = 208686;

	double yMax = 1.5, yMin = 0.0;

	TLine blinedown(bdown-2,yMin,bdown-2,yMax);
	TLine blineup(bup-2,yMin,bup-2,yMax);

	TLine clinedown(cdown-2,yMin,cdown-2,yMax);
	TLine clineup(cup-2,yMin,cup-2,yMax);

	TLine dlinedown(ddown-2,yMin,ddown-2,yMax);
	TLine dlineup(dup-2,yMin,dup-2,yMax);

	blinedown.SetLineColor(kRed);
	blinedown.SetLineWidth(1);

	blineup.SetLineColor(kRed);
	blineup.SetLineWidth(1);

	clinedown.SetLineColor(kRed);
	clinedown.SetLineWidth(1);

	clineup.SetLineColor(kRed);
	clineup.SetLineWidth(1);

	dlinedown.SetLineColor(kRed);
	dlinedown.SetLineWidth(1);

	dlineup.SetLineColor(kRed);
	dlineup.SetLineWidth(1);

	TH1F* hlt4 = (TH1F*)f->Get("lumiJPsirate_HLT4");
	TH1F* hlt8 = (TH1F*)f->Get("lumiJPsirate_HLT8");

	TH1F* lumihlt4 = (TH1F*)f->Get("JPsilumi_rebin_HLT4");
	TH1F* lumihlt8 = (TH1F*)f->Get("JPsilumi_rebin_HLT8");

	TH1F* all_runs_hlt4 = new TH1F("allruns_hlt4", "JPsis from B^{0}_{s} candidate; Run Number;no. of J/Psi per lumi (mb)",18000, 193000, 210000);
	TH1F* all_runs_hlt8 = new TH1F("allruns_hlt8", "allruns_hlt8; Run[#];no. of J/Psi per lumi (mb)",18000, 193000, 210000);

	TH1F* b_runs_hlt4 = new TH1F("b_runs_hlt4", "JPsis from B^{0}_{s} candidate - Dataset B; Run Number;no. of J/Psi per lumi (mb)",4500, 193500, 197000);
	TH1F* b_runs_hlt8 = new TH1F("b_runs_hlt8", "b_runs_hlt8; Run[#];J/Psi[#/mb]",4500, 193500, 197000);

	TH1F* c_runs_hlt4 = new TH1F("c_runs_hlt4", "JPsis from B^{0}_{s} candidate - Dataset C; Run Number;no. of J/Psi per lumi (mb)",8000, 197000, 205000);
	TH1F* c_runs_hlt8 = new TH1F("c_runs_hlt8", "c_runs_hlt8; Run[#];J/Psi[#/mb]",8000, 197000, 205000);

	TH1F* d_runs_hlt4 = new TH1F("d_runs_hlt4", "JPsis from B^{0}_{s} candidate - Dataset D; Run Number;no. of J/Psi per lumi (mb)",8000, 202000, 210000);
	TH1F* d_runs_hlt8 = new TH1F("d_runs_hlt8", "d_runs_hlt8; Run[#];J/Psi[#/mb]",8000, 202000, 210000);

	TCanvas canvas("canvas","canvas",1200,800);

	TLegend legend(0.6,0.8,0.9,0.9);

	b_runs_hlt4->SetLineColor(1);
	b_runs_hlt4->SetLineStyle(0);
	b_runs_hlt4->SetLineWidth(1);
	b_runs_hlt4->SetMarkerStyle(20);
	b_runs_hlt4->SetMarkerColor(kRed);
	b_runs_hlt4->SetMaximum(yMax);
	b_runs_hlt4->SetMinimum(yMin);

	b_runs_hlt8->SetLineColor(1);
	b_runs_hlt8->SetLineStyle(0);
	b_runs_hlt8->SetLineWidth(1);
	b_runs_hlt8->SetMarkerStyle(20);
	b_runs_hlt8->SetMarkerColor(kBlue);
	b_runs_hlt8->SetMaximum(yMax);
	b_runs_hlt8->SetMinimum(yMin);

	c_runs_hlt4->SetLineColor(1);
	c_runs_hlt4->SetLineStyle(0);
	c_runs_hlt4->SetLineWidth(1);
	c_runs_hlt4->SetMarkerStyle(20);
	c_runs_hlt4->SetMarkerColor(kRed);
	c_runs_hlt4->SetMaximum(yMax);
	c_runs_hlt4->SetMinimum(yMin);

	c_runs_hlt8->SetLineColor(1);
	c_runs_hlt8->SetLineStyle(0);
	c_runs_hlt8->SetLineWidth(1);
	c_runs_hlt8->SetMarkerStyle(20);
	c_runs_hlt8->SetMarkerColor(kBlue);
	c_runs_hlt8->SetMaximum(yMax);
	c_runs_hlt8->SetMinimum(yMin);


	d_runs_hlt4->SetLineColor(1);
	d_runs_hlt4->SetLineStyle(0);
	d_runs_hlt4->SetLineWidth(1);
	d_runs_hlt4->SetMarkerStyle(20);
	d_runs_hlt4->SetMarkerColor(kRed);
	d_runs_hlt4->SetMaximum(yMax);
	d_runs_hlt4->SetMinimum(yMin);

	d_runs_hlt8->SetLineColor(1);
	d_runs_hlt8->SetLineStyle(0);
	d_runs_hlt8->SetLineWidth(1);
	d_runs_hlt8->SetMarkerStyle(20);
	d_runs_hlt8->SetMarkerColor(kBlue);
	d_runs_hlt8->SetMaximum(yMax);
	d_runs_hlt8->SetMinimum(yMin);

	all_runs_hlt4->SetLineColor(1);
	all_runs_hlt4->SetLineStyle(0);
	all_runs_hlt4->SetLineWidth(1);
	all_runs_hlt4->SetMarkerStyle(20);
	all_runs_hlt4->SetMarkerColor(kRed);
	all_runs_hlt4->SetMaximum(yMax);
	all_runs_hlt4->SetMinimum(yMin);

	all_runs_hlt8->SetLineColor(1);
	all_runs_hlt8->SetLineStyle(0);
	all_runs_hlt8->SetLineWidth(1);
	all_runs_hlt8->SetMarkerStyle(20);
	all_runs_hlt8->SetMarkerColor(kBlue);
	all_runs_hlt8->SetMaximum(yMax);
	all_runs_hlt8->SetMinimum(yMin);


	for(int i = 0; i<hlt4->GetNbinsX();++i)
		{
			double center = hlt4->GetBinCenter(i);
			std::cout << center << std::endl;
			if(center >= bdown && center<=bup)
				{
					b_runs_hlt4->SetBinContent(b_runs_hlt4->FindBin(center),hlt4->GetBinContent(i));
					b_runs_hlt4->SetBinError(b_runs_hlt4->FindBin(center),hlt4->GetBinError(i));
				}
			else if(center >= cdown && center<=cup)
				{
					c_runs_hlt4->SetBinContent(c_runs_hlt4->FindBin(center),hlt4->GetBinContent(i));
					c_runs_hlt4->SetBinError(c_runs_hlt4->FindBin(center),hlt4->GetBinError(i));
				}
			else if(center >= ddown && center<=dup)
				{
					d_runs_hlt4->SetBinContent(d_runs_hlt4->FindBin(center),hlt4->GetBinContent(i));
					d_runs_hlt4->SetBinError(d_runs_hlt4->FindBin(center),hlt4->GetBinError(i));
				}

				all_runs_hlt4->SetBinContent(all_runs_hlt4->FindBin(center),hlt4->GetBinContent(i));
				all_runs_hlt4->SetBinError(all_runs_hlt4->FindBin(center),hlt4->GetBinError(i));

		}

		for(int i = 0; i<hlt8->GetNbinsX();++i)
			{
				double center = hlt8->GetBinCenter(i);
				std::cout << center << std::endl;
				if(center >= 193833 && center<=196531)
					{
						b_runs_hlt8->SetBinContent(b_runs_hlt8->FindBin(center),hlt8->GetBinContent(i));
						b_runs_hlt8->SetBinError(b_runs_hlt8->FindBin(center),hlt8->GetBinError(i));
					}
					else if(center >= cdown && center<=cup)
						{
							c_runs_hlt8->SetBinContent(c_runs_hlt8->FindBin(center),hlt8->GetBinContent(i));
							c_runs_hlt8->SetBinError(c_runs_hlt8->FindBin(center),hlt8->GetBinError(i));
						}
					else if(center >= ddown && center<=dup)
						{
							d_runs_hlt8->SetBinContent(d_runs_hlt8->FindBin(center),hlt8->GetBinContent(i));
							d_runs_hlt8->SetBinError(d_runs_hlt8->FindBin(center),hlt8->GetBinError(i));
						}

						all_runs_hlt8->SetBinContent(all_runs_hlt8->FindBin(center),hlt8->GetBinContent(i));
						all_runs_hlt8->SetBinError(all_runs_hlt8->FindBin(center),hlt8->GetBinError(i));

			}

		b_runs_hlt4->Draw("PE");
		b_runs_hlt8->Draw("PEsame");

		legend.AddEntry("b_runs_hlt4","HLT_4 + pt cut","lep");
		legend.AddEntry("b_runs_hlt8","HLT_8 + Lxy cut","lep");
		legend.Draw();

		canvas.SaveAs("bRun.eps");

		c_runs_hlt4->Draw("PE");
		c_runs_hlt8->Draw("PEsame");
		legend.Draw();
		canvas.SaveAs("cRun.eps");

		d_runs_hlt4->Draw("PE");
		d_runs_hlt8->Draw("PEsame");
		legend.Draw();
		canvas.SaveAs("dRuns.eps");

		all_runs_hlt4->Draw("PE");
		//all_runs_hlt4->SetOptStat(0000);
		all_runs_hlt8->Draw("PEsame");

		blinedown.Draw();
		blineup.Draw();

		dlinedown.Draw();
		dlineup.Draw();

		clinedown.Draw();
		legend.Draw();
		// clineup.Draw();

		canvas.SaveAs("allRuns.eps");

}
