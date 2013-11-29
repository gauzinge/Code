//  
//
//  Created by Georg Auzinger on 15.03.13.
//
//
#include <math.h>
#include <stdio.h>
#include <vector>
#include <set>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>

// #include <dirent.h>
//ROOT
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TFitResultPtr.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TFitResult.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMatrixT.h>

const double PI  =3.141592653589793238462;

struct event{
    TFile* hitfile;
    TTree* hittree;
	int _EvtNr;
	int _RunNr;
	
	double _x[6];
	double _y[6];
	double _z[6];
    void readtree(std::string);
    int get_n_entries();
    void get_entry(int);
};

void event::readtree(std::string filename)
{
	//first read the tree
	this->hitfile = TFile::Open(filename.c_str(),"OPEN");
	if (!hitfile) std::cerr << "Could not open Hit File!" << std::endl;
	else
	{
		hitfile->cd("MyEUTelFitTuple");
		this->hittree = (TTree*)gDirectory->Get("EUFit");
		if (!hittree) std::cerr << "Could not open Hit Tree!" << std::endl;
		else 
		{
			TBranch* eventbranch = 0;
			hittree->SetBranchAddress("EvtNr",&_EvtNr,&eventbranch);
			
			TBranch* runbranch = 0;
			hittree->SetBranchAddress("RunNr",&_RunNr,&runbranch);
			
			TBranch* x_0_branch = 0;
			TBranch* y_0_branch = 0;
			TBranch* z_0_branch = 0;
			hittree->SetBranchAddress("measX_0",&_x[0],&x_0_branch);
			hittree->SetBranchAddress("measY_0",&_y[0],&y_0_branch);
			hittree->SetBranchAddress("measZ_0",&_z[0],&z_0_branch);
			
			TBranch* x_1_branch = 0;
			TBranch* y_1_branch = 0;
			TBranch* z_1_branch = 0;
			hittree->SetBranchAddress("measX_1",&_x[1],&x_1_branch);
			hittree->SetBranchAddress("measY_1",&_y[1],&y_1_branch);
			hittree->SetBranchAddress("measZ_1",&_z[1],&z_1_branch);
			
			TBranch* x_2_branch = 0;
			TBranch* y_2_branch = 0;
			TBranch* z_2_branch = 0;
			hittree->SetBranchAddress("measX_2",&_x[2],&x_2_branch);
			hittree->SetBranchAddress("measY_2",&_y[2],&y_2_branch);
			hittree->SetBranchAddress("measZ_2",&_z[2],&z_2_branch);
			
			TBranch* x_3_branch = 0;
			TBranch* y_3_branch = 0;
			TBranch* z_3_branch = 0;
			hittree->SetBranchAddress("measX_3",&_x[3],&x_3_branch);
			hittree->SetBranchAddress("measY_3",&_y[3],&y_3_branch);
			hittree->SetBranchAddress("measZ_3",&_z[3],&z_3_branch);
			
			TBranch* x_4_branch = 0;
			TBranch* y_4_branch = 0;
			TBranch* z_4_branch = 0;
			hittree->SetBranchAddress("measX_4",&_x[4],&x_4_branch);
			hittree->SetBranchAddress("measY_4",&_y[4],&y_4_branch);
			hittree->SetBranchAddress("measZ_4",&_z[4],&z_4_branch);
			
			TBranch* x_5_branch = 0;
			TBranch* y_5_branch = 0;
			TBranch* z_5_branch = 0;
			hittree->SetBranchAddress("measX_5",&_x[5],&x_5_branch);
			hittree->SetBranchAddress("measY_5",&_y[5],&y_5_branch);
			hittree->SetBranchAddress("measZ_5",&_z[5],&z_5_branch);
		}
	}
}

int event::get_n_entries() {
    return static_cast<int>(hittree->GetEntries());
}

void event::get_entry(int i_entry){
    hittree->GetEntry(i_entry);
}

struct track{
	int EvtNr;
	int RunNr;

	// Track Parameters a+b*z
	double aX;
	double aY;
	double bX;
	double bY;

	// Track Parameter Errors propagate as sigma(z)=sqrt(sigma_a^2 + sigma_b^2*z^2 + 2*sigma_ab) for both coordinates yields the track uncertanty sigma_x,y(z)
	double sigma_aX;
	double sigma_aY;
	double sigma_bX;
	double sigma_bY;
};

double get_fitrange(TH1D* h) {
    double maxval = h->GetMaximum();
    double firstbin = h->GetBinCenter(h->FindFirstBinAbove(maxval/2));
    double lastbin = h->GetBinCenter(h->FindLastBinAbove(maxval/2));
    double fwhm = (lastbin-firstbin);
    return (2*fwhm/2);
}


void get_track_uncert(std::string filename, double *errorsX, double *errorsY)
{	
	//instantiate event
	event myevent;
	myevent.readtree(filename);
	
	// her should be the loop for the iteration 
	//now exclude 1 plane at a time
	int nbins=2000;
	double minbin=-.15;
	double maxbin=.15;

	std::vector<TH2D*> sigmahisto;
	std::vector<TH2D*> residualhisto;
	std::vector<TH2D*> sigmaDUT;
	
	//prepare histogram vectors
	for (int excluded=0 ; excluded<6 ; excluded++)
	{
		//std::cout the used errors for this iteration:
		std::cout << "Input errors plane " << excluded << " x " << errorsX[excluded] << " y " << errorsY[excluded] << std::endl;
				
		std::stringstream ss;
		ss << excluded;
		std::string sigmaHistoName="Sigma_plane_" + ss.str();
		std::string residualHistoName="Residuals_plane_"+ss.str();
		std::string sigmaDUTName="Sigma_DUT_plane_" + ss.str();
		
		TH2D* dummy1 = new TH2D(sigmaHistoName.c_str(),sigmaHistoName.c_str(),nbins,minbin,maxbin,nbins,minbin,maxbin);
		
		TH2D* dummy2 = new TH2D(residualHistoName.c_str(),residualHistoName.c_str(),nbins,minbin,maxbin,nbins,minbin,maxbin);
		
		TH2D* dummy3 = new TH2D(sigmaDUTName.c_str(),sigmaDUTName.c_str(),nbins,-0.1,0.1,nbins,-0.1,0.1);
		
		sigmahisto.push_back(dummy1);
		residualhisto.push_back(dummy2);
		sigmaDUT.push_back(dummy3);
	}
	
	//finally, DATA
	for (int excluded=0 ; excluded<6 ; excluded++)
	{
		std::cout << "Excluding Plane " << excluded << std::endl;
		
		for (int events=0; events<myevent.get_n_entries(); events++)
		{
			myevent.get_entry(events);
						
			double residualX=0;
			double residualY=0;
			//here i have all hits per event
			
			//TGraphErrors for Fit!
			TGraphErrors* grx = new TGraphErrors();
			TGraphErrors* gry = new TGraphErrors();
			
			//fill them with hits from all planes but the excluded
			for (int plane=0 ; plane<6 ; plane++)
			{
				if (plane!=excluded) //this is a valid plane, eventually sanity check
				{
					int npx = grx->GetN();
					int npy = gry->GetN();
					
					grx->SetPoint(npx,myevent._z[plane],myevent._x[plane]);
					grx->SetPointError(npx,0,errorsX[plane]);
					
					gry->SetPoint(npy,myevent._z[plane],myevent._y[plane]);
					gry->SetPointError(npy,0,errorsY[plane]);
				}
				else
				{
					residualX=myevent._x[plane];
					residualY=myevent._y[plane];
				}
			}
			//graphs for each event filled with data from 5 planes
			
			TFitResultPtr resx = grx->Fit("pol1","SQ");
			double x = resx.Get()->Parameter(0)+resx.Get()->Parameter(1)*myevent._z[excluded];
			double sigma_x = TMath::Sqrt(resx.Get()->ParError(0)*resx.Get()->ParError(0)+resx.Get()->ParError(1)*resx.Get()->ParError(1)*myevent._z[excluded]*myevent._z[excluded]+2*(resx.Get()->GetCovarianceMatrix()(0,1)*myevent._z[excluded]));
				
			TFitResultPtr resy = gry->Fit("pol1","SQ");
			double y = resy.Get()->Parameter(0)+resy.Get()->Parameter(1)*myevent._z[excluded];
			double sigma_y = TMath::Sqrt(resy.Get()->ParError(0)*resy.Get()->ParError(0)+resy.Get()->ParError(1)*resy.Get()->ParError(1)*myevent._z[excluded]*myevent._z[excluded]+2*(resy.Get()->GetCovarianceMatrix()(0,1)*myevent._z[excluded]));
			 
			//now Fill 2D histograms
			sigmahisto.at(excluded)->Fill(sigma_x,sigma_y); //should be same
			residualhisto.at(excluded)->Fill(residualX-x,residualY-y); //residual distribution
		} //end of event loop
	} //end of exclusion loop
			
	double newerrorsX[6]={0,0,0,0,0,0};
	double newerrorsY[6]={0,0,0,0,0,0};
			
	//Here I could Plot the residuals, why not actually
	TCanvas* rescanvas = new TCanvas("residuals","residuals");
	rescanvas->Divide(3,2);
	
	for (int i=0 ; i<6 ; i++)
	{
		rescanvas->cd(i+1);
		residualhisto.at(i)->Draw();
		
		//for every plane, compute the intrinsic resolution which is 
		//sqrt(sigma_residuals^2-sigma_track^2)
		
		newerrorsX[i] = residualhisto.at(i)->GetRMS(1)/(TMath::Sqrt(1+((sigmahisto.at(i)->GetMean(1)*sigmahisto.at(i)->GetMean(1))/(errorsX[i]*errorsX[i]))));
		newerrorsY[i] = residualhisto.at(i)->GetRMS(2)/(TMath::Sqrt(1+(sigmahisto.at(i)->GetMean(2)*sigmahisto.at(i)->GetMean(2)/errorsY[i]*errorsY[i])));
		
		//output to debug
		std::cout << "Intrinsic Sensor Resolution [mum] for plane " << i << " in x " << newerrorsX[i] << " and y " << newerrorsY[i] << std::endl;
		
		errorsX[i]=newerrorsX[i];
		errorsY[i]=newerrorsY[i];
	}
}

void fit_hits(event myevent, track &mytrack, double *errorsX, double *errorsY, std::set<int> excluded_planes)
{
	//TGraphErrors for Fit!
	TGraphErrors* grx = new TGraphErrors();
	TGraphErrors* gry = new TGraphErrors();

	//fill them with hits from all planes
	for (int plane=0 ; plane<6 ; plane++)
	{
		if (excluded_planes.find(plane) == excluded_planes.end())
		{
			int npx = grx->GetN();
			int npy = gry->GetN();
	
			grx->SetPoint(npx,myevent._z[plane],myevent._x[plane]);
			grx->SetPointError(npx,0,errorsX[plane]);
	
			gry->SetPoint(npy,myevent._z[plane],myevent._y[plane]);
			gry->SetPointError(npy,0,errorsY[plane]);
		}
	}

	//Fit
	TFitResultPtr resx = grx->Fit("pol1","SQ");
	TFitResultPtr resy = gry->Fit("pol1","SQ");

	//set members
	mytrack.EvtNr = myevent._EvtNr;
	mytrack.RunNr = myevent._RunNr;
	
	// Track Parameters a+b*z
	mytrack.aX = resx.Get()->Parameter(0);
	mytrack.bX = resx.Get()->Parameter(1);
	
	mytrack.sigma_aX = resx.Get()->ParError(0);
	mytrack.sigma_bX = resx.Get()->ParError(1);

	mytrack.aY = resy.Get()->Parameter(0);
	mytrack.bY = resy.Get()->Parameter(1);
	
	mytrack.sigma_aY = resy.Get()->ParError(0);
	mytrack.sigma_bY = resy.Get()->ParError(1);
}


void write_tracks(std::string filename, double *errorsX, double *errorsY, std::set<int> excluded_planes)
{	
	track mytrack;
	//create track file and tree
	std::string trackfilename = filename.substr(0,9) + "_tracks.root";
	TFile* trackfile = TFile::Open(trackfilename.c_str(),"RECREATE");
	if (!trackfile) std::cerr << "Could not open Trackfile!" << std::endl;
	else
	{
		TTree* tracktree = new TTree("TrackParameters","TrackParameters");
		if (!tracktree) std::cerr << "Could not create Track Tree!" << std::endl;
		else
		{
			//define branches and point them to members of mytrack
			tracktree->Branch("EvtNr",&mytrack.EvtNr);
			tracktree->Branch("RunNr",&mytrack.RunNr);
			tracktree->Branch("aX",&mytrack.aX); //intercept x
			tracktree->Branch("bX",&mytrack.bX); //slope x
			tracktree->Branch("sigma_aX",&mytrack.sigma_aX); //error ax
			tracktree->Branch("sigma_bX",&mytrack.sigma_bX); //error bx
			tracktree->Branch("aY",&mytrack.aY); //intercept y
			tracktree->Branch("bY",&mytrack.bY); //slope y
			tracktree->Branch("sigma_aY",&mytrack.sigma_aY); //error ay
			tracktree->Branch("sigma_bY",&mytrack.sigma_bY); //error by
	
			//instantiate event
			event myevent;
			myevent.readtree(filename);
	
			//loop over events
			for (int events=0; events<myevent.get_n_entries(); events++)
			{
				myevent.get_entry(events);
				fit_hits(myevent, mytrack, errorsX, errorsY, excluded_planes);
				tracktree->Fill();
			}
			myevent.hitfile->Close();
		}
		trackfile->cd();
		tracktree->Write();
		trackfile->Write();
		trackfile->Close();
	}
}

void track_converter(int runnumber, int exclude=5, int nsteps=3)
{
	std::stringstream filestream;
	filestream << "run" << setfill('0') << setw(6) << runnumber << "-fitter.root";
	std::string filename = filestream.str();
	
	std::set<int> excluded_planes;
	excluded_planes.insert(exclude);
	
	//initial measurement errors, will be overwritten in any iteration
	double errorsX[6]={0.01,0.01,0.01,0.01,0.01,0.01};
	double errorsY[6]={0.01,0.01,0.01,0.01,0.01,0.01};
	
	std::cout << "Computing Resolution of individual Telescope Planes in " << nsteps << " iterations!" << std::endl;
	
	//initialize std::vector of TGraph* for visualization of convergence
	std::vector<TGraph*> xvector;
	std::vector<TGraph*> yvector;
	
	for (int plane=0 ; plane<6 ; plane++)
	{
		std::stringstream xstream;
		std::stringstream ystream;
		
		xstream << "Resolution plane " << plane;
		ystream << "Resolution plane " << plane;
		
		TGraph* dummygraphx = new TGraph();
		TGraph* dummygraphy = new TGraph();
		
		dummygraphx->SetTitle(xstream.str().c_str());
		dummygraphy->SetTitle(ystream.str().c_str());
		
		dummygraphx->SetLineColor(1);
		dummygraphy->SetLineColor(2);
		
		dummygraphx->GetYaxis()->SetRangeUser(0.001,0.05);
		dummygraphy->GetYaxis()->SetRangeUser(0.001,0.05);
		
		xvector.push_back(dummygraphx);
		yvector.push_back(dummygraphy);
	}
	
	for (int iterations=0 ; iterations<nsteps ; iterations++)
	{
		std::cout << "Running iteration " << iterations << std::endl;
		get_track_uncert(filename, errorsX, errorsY);
		
		if (iterations != 0)
		{
			for (int plane=0 ; plane<6; plane++)
			{
				xvector.at(plane)->SetPoint(xvector.at(plane)->GetN(),iterations,errorsX[plane]);
				yvector.at(plane)->SetPoint(yvector.at(plane)->GetN(),iterations,errorsY[plane]);
			}
		}	
	}
	
	TCanvas* iterationcanvas = new TCanvas("iterationcanvas","iterationcanvas");
	iterationcanvas->Divide(3,2);
	for (int plane=0 ; plane<6 ; plane++)
	{
		iterationcanvas->cd(plane+1);
		TMultiGraph* multigraph = new TMultiGraph("","Resolution vs. Iteration #;Iteration #;Resolution [mm]");
		multigraph->Add(xvector.at(plane),"PL");
		multigraph->Add(yvector.at(plane),"PL");
		// xvector.at(plane)->Draw("APL");
// 		yvector.at(plane)->Draw("PL same");
		multigraph->Draw("A");
	}
	
	//now i have the individual sensor plane resolutions (that are the error of the measurement) to assign to the measured hits -> TGraphErrors -> Fit pol1 -> extract parameters and parameter errors -> write to tree
	
	std::cout << "Parametrizing Tracks with planes ";
	for (std::set<int>::iterator it = excluded_planes.begin() ; it!=excluded_planes.end() ; it++) std::cout << *it << " ";
	std::cout << " excluded!" << std::endl;
	
	write_tracks(filename, errorsX, errorsY, excluded_planes);
}