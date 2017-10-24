#include <TH2F.h>
#include <TFile.h>
#include <iostream>
#include <map>
#include <vector>
#include <algorithm>


void jpsiAnalize(TFile* f)
{
	TH2F* h = (TH2F*)f->Get("runLumi");
	std::map<int, std::vector < int > > runs;	
	std::map<int, std::vector < int > >::iterator runsIt;
	int nBinY = h->GetNbinsY();
	int nBinX = h->GetNbinsX();
		
	std::vector <int> lumis;
	double x, y;
	//std::cout << "{";
	double prevY=0;
	double lumicount=0;
	double lumMax = 0;
	for(int i = 0; i < nBinX + 4 ;++i)
	{
		prevY=0.0;
		double runcounter = 0;
		//x = h->GetXaxis()->GetBinCenter(i);
		//std::cout << "\""<<x<<"\":[";
		double candcounter=0.;
		for(int j = 0; j< nBinY + 4; ++j)
			{
			
				int bin = h->GetBin(i,j);
				double cont = h->GetBinContent(bin);
			
				x = h->GetXaxis()->GetBinCenter(i);
				candcounter +=cont;
			
				if(cont>0.0)
				{	
					runcounter++;					
					runs[int(x)] = lumis; 
					y = h->GetYaxis()->GetBinCenter(j);	
					if(prevY<y-1 && prevY!=0.0)
					std::cout << x << " " << y << " " << prevY <<std::endl;
					prevY=y;
				}
							
			}
		//std::cout<<"], ";
		if(runcounter >0.0)		
		std::cout << x <<" -> run counter " << runcounter << " -> cand counter :" << candcounter <<std::endl; 
		lumMax = std::max(lumMax,runcounter);
	}
	//std::cout<<"}";
	
	for(int i = 0; i < nBinX + 4 ;++i)
	{
		 for(int j = 0; j< nBinY + 4; ++j)
		 {
			 x = h->GetXaxis()->GetBinCenter(i) ;
			 int bin = h->GetBin(i,j);
			 double cont = h->GetBinContent(bin);
		 	 if(cont>0.0)
			 {
				 y = h->GetYaxis()->GetBinCenter(j);  
				 runs[int(x)].push_back(int(y));
				 ++lumicount;
			 }

		 }
	}

	std::cout << "{";	
	for(runsIt=runs.begin();runsIt!=runs.end();runsIt++)
	{

		std::cout << "\""<<runsIt->first<<"\":[";
		int prevLum=0;
		for(int k = 0; k<runsIt->second.size();k++)
		{
			if(k==0)
			{
				std::cout<<"[";
				std::cout<<(runsIt->second)[k]<<", ";	
			
				if(k==runsIt->second.size()-1)
			                        {
					                                std::cout<<(runsIt->second)[k]<<"]";
							                        }

			}
			else
			if(!((runsIt->second)[k]==prevLum+1))
			{
				std::cout<<prevLum<<"], ["<<(runsIt->second)[k]<<", ";		
				if(k==runsIt->second.size()-1)
				{
					std::cout<<(runsIt->second)[k]<<"]";
				}
				
			}
			else
			{
				if(k==runsIt->second.size()-1)
					                                {
										                                        std::cout<<(runsIt->second)[k]<<"]";
															                                }
			}

			

			prevLum = (runsIt->second)[k];
		}	
			
		
		std::cout << "], ";
	}
	std::cout<<"}";
	std::cout<<std::endl;
	std::cout<<lumicount<<std::endl;
	std::cout<<lumMax<<std::endl;
}

