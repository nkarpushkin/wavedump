#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <array>
#include <functional>
#include <iostream>
//#include <string_view>

#include <TFile.h>
#include <TNtuple.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TTree.h>

#include <TCanvas.h>
#include <TH1F.h>

#include <TLine.h>

#include "../PronyFitter.cxx"
#include "../PronyFitter.h"
using namespace std;

#define EVENT_SIZE 1024
#define N_CHANNELS 16

int fit_channel = 6;
int com_channel1 = 4;
int com_channel2 = 8;

typedef struct _CaenEvent
{
  int recLength;      // Record Length: 1024
  int boardID;        // BoardID: 31
  int channel;        // Channel: 0
  int eventNumber;    // Event Number: 0
  int pattern;        // Pattern: 0x0000
  int timeStamp;      // Trigger Time Stamp: 392909414
  int dcOffset;       // DC offset (DAC): 0x7FFF
  int startIndexCell; // Start Index Cell: 562
  short int data[1024];
  float dataInt[1024];
  float dataIntegral;
  short int max;
  short int min;
  short int minPed;
  float avgPed;
  short int avg;
  float timeAtThresh;
  float timeOverThresh;
} CaenEvent;

enum _EventHeaderFields
{
  REC_LENGTH = 0,
  BOARD_ID,
  CHANNEL,
  EVENT_NUMBER,
  PATTERN,
  TIME_STAMP,
  DC_OFFSET,
  START_INDEX_CELL,

  EVENT_FIELDS_NUMBER
} EventFields;

const char EventHeaderFieldStrings[][100] = {
    "Record Length",
    "BoardID",
    "Channel",
    "Event Number",
    "Pattern",
    "Trigger Time Stamp",
    "DC offset (DAC)",
    "Start Index Cell"};

typedef std::vector<CaenEvent> VectorOfCaenEvents;
typedef std::vector<VectorOfCaenEvents> VectorOfVectorOfCaenEvents;

void MeanRMScalc(short *DataArr, float* Mean, float* RMS, int begin, int end, int step = 1)
{
    begin = (begin < 0)? 0 : begin;

    if(begin > end){float swap=end;end=begin;begin=swap;};
    step = TMath::Abs(step);

    *Mean = *RMS = 0.; int Delta = 0;
    for(int n=begin; n<=end; n += step){ *Mean += DataArr[n]; Delta++;}
    *Mean /= (float)Delta;

    for(int n=begin; n<=end; n += step) *RMS += (DataArr[n] - *Mean) * (DataArr[n] - *Mean);
    *RMS = TMath::Sqrt( *RMS/((float)Delta) );

    //printf("AMPL %.2f, RMS %.2f\n",*Mean,*RMS);
}

float LevelBy2Points(float X1, float Y1, float X2, float Y2, float Y0)
{
    //            [X2, Y2] 0
    //                   *
    //                 *
    //   Y0--------- 0
    //             *  X0 (returned)
    //  [X1, Y1] 0
    return (X1*Y0 - X1*Y2 - X2*Y0 + X2*Y1) / (Y1-Y2);
}
float GoToLevel(short *DataArr, float Level, int *point, int iterator, int iLastPoint)
{

    float ResultTime;
    while( (*point>=0)&&(*point<iLastPoint) )
    {

	if( (Level-DataArr[*point])*(Level-DataArr[*point+iterator]) <= 0  ){
	    ResultTime = LevelBy2Points((float)(*point),DataArr[*point],(float)(*point+iterator),DataArr[*point+iterator],Level);
	    return ResultTime;
	}
	//printf("point %i, ampl1 %.2f ampl2 %.2f\n",*point,Level-DataArr[*point],Level-DataArr[*point+iterator]);
	*point += iterator;
    } //while

    *point = -1;
    return 0;
}

void saveHistogramsToRootFile(std::vector<TCanvas *> &canvases, const std::string &filename)
{
  // Create a ROOT file to store the histograms
  TFile file(filename.c_str(), "RECREATE");

  // Loop over the canvases and save their histograms
  for (auto canvas : canvases)
  {
    canvas->Write();

    for (int j = 0; j < canvas->GetListOfPrimitives()->GetSize(); j++)
    {
      TObject *obj = canvas->GetListOfPrimitives()->At(j);
      if (obj->InheritsFrom("TPad"))
      {
        TPad *pad = (TPad *)obj;
        for (int k = 0; k < pad->GetListOfPrimitives()->GetSize(); k++)
        {
          TObject *subobj = pad->GetListOfPrimitives()->At(k);
          if (subobj->InheritsFrom("TH1"))
          {
            TH1 *hist = (TH1 *)subobj;
            hist->Write();
          }
        }
      }
    }
  }

  file.Close();
}

void wavedump2024(char *folder = 0, int nEventsToProcess = -1,
                              char *outFileName = "data.root",
                              int signalBegin = 400, int signalEnd = 900,
                              int pedestalBegin = 10, int pedestalEnd = 210,
                              float constantFractionThresh = 0.2,
                              float relThresh = 3000., int showData = 0, int eventToShow = -1,
                              const int N_CHANNELS_TO_SHOW = 4,
                              Double_t rangeUser1 = 0, Double_t rangeUser2 = 4000)
{

  TTree *_bigTree = new TTree("bigtree", "bigtree");
  TString pics_file_name = "picresults_wavedump/" + TString(outFileName);
  std::vector<TCanvas *> canvas_array;

  int _treeFill_iEvt;
  _bigTree->Branch("iEvt", &_treeFill_iEvt, "iEvt/I");

  float _treeFill_adc[16];
  float _treeFill_qdc[16];
  float _treeFill_tdc[16];
  float _treeFill_tot[16];
  float _treeFill_ped[16];
  _bigTree->Branch("adc", _treeFill_adc, "adc[16]/F");
  _bigTree->Branch("qdc", _treeFill_qdc, "qdc[16]/F");
  _bigTree->Branch("tdc", _treeFill_tdc, "tdc[16]/F");
  _bigTree->Branch("tot", _treeFill_tot, "tot[16]/F");
  _bigTree->Branch("ped", _treeFill_ped, "ped[16]/F");

  TCanvas *c1 = 0;
  TH1F *hf1 = 0;

  //##################################
  const int model_order = 5;
  const int exponents = 2;
  const float fit_R2_thr = 0.02;
  //##################################
  int counter = 0;
  std::vector<TH2F*> harmHistVect; 
  harmHistVect.resize(exponents);
  for(int i =0; i < exponents; i++)
      harmHistVect.at(i) = new TH2F(Form("harmonic%i", i+1),Form("harmonic%i", i+1), 500,0,1.5,500,-1,1);

  std::vector<TH1F*> tauHistVect; 
  tauHistVect.resize(exponents);
  for(int i =0; i < exponents; i++)
      tauHistVect.at(i) = new TH1F(Form("tau%i", i+1),Form("tau%i; tau [ns]; Counts[]", i+1), 500,0,20);

  std::vector<TH2F*> amplCorrVect; 
  amplCorrVect.resize(exponents*exponents);
  int c_corr = 0;
  for(int i =0; i < exponents; i++)
    for(int j =0; j < exponents; j++){
      amplCorrVect.at(c_corr) = new TH2F(Form("ampl%i_vs_ampl%i", i+1, j+1),Form("ampl%i_vs_ampl%i", i+1, j+1), 1500,-3000,3000,1500,-3000,3000);
      c_corr++;
    }

  TH1D* wfm_ethalon = new TH1D("wfm_ethalon", "wfm_ethalon", 1024,0,1024);

 // FIAN
 
  int signalType[16] = {-1, -1, -1, -1,
                        -1, 1, 1, 1,
                         1, -1, 1, -1,
                        -1, -1, -1, -1};
  
/*
  // two scintillators cosmic. "Two channel" readout
  int signalType[16] = {1, 1, -1, -1,
                        -1, 1, 1, 1,
                        -1, -1, -1, -1,
                        -1, -1, -1, -1};
*/

  char path[1000];
  std::ifstream datafiles[N_CHANNELS];
  for (int i = 0; i < N_CHANNELS; i++)
  {
    if (folder == 0)
    {
      sprintf(path, "./wave_%d.dat", i);
      clog << "Opening file: " << path << endl;
      datafiles[i].open(path);
    }
    else
    {
      sprintf(path, "%s/wave_%d.dat", folder, i);
      clog << "Opening file: " << path << endl;
      datafiles[i].open(path);
    }
  }

  if (datafiles[0].is_open())
  {

    TFile outfile(outFileName, "RECREATE", "data", 5);
    std::vector<std::vector<CaenEvent>> data; // = new std::vector<std::vector<CaenEvent> >();

    datafiles[0].seekg(0, datafiles[0].end);
    long int fileLength = datafiles[0].tellg();
    datafiles[0].seekg(0, datafiles[0].beg);

    cout << "Number of events = " << fileLength / (4 * 1024) << endl;

    long int nEvents = fileLength / (4 * 1024);

    if (nEventsToProcess > 0 && nEventsToProcess <= nEvents)
      nEvents = nEventsToProcess;

    cout << "Number of events to process= " << nEvents << endl;

    for (int iEvt = 0; iEvt < nEvents; iEvt++)
    {

      _treeFill_iEvt = iEvt;
      for (int k = 0; k < 16; k++)
      {
        _treeFill_adc[k] = 0;
        _treeFill_qdc[k] = 0;
        _treeFill_tdc[k] = 0;
        _treeFill_tot[k] = 0;
      }

      std::vector<CaenEvent> channelData;

      for (int iCh = 0; iCh < N_CHANNELS; iCh++)
      {

        CaenEvent ev;

        ev.recLength = 1024;
        ev.eventNumber = iEvt;

        ev.min = SHRT_MAX; // no typo
        ev.max = SHRT_MIN; // no typo
        ev.minPed = SHRT_MAX; // no typo

        int avg = 0;
        float avgPed = 0;
        float buf;
        float dataInt = 0; // for integral
        float q = 0;

        int pedBegin = pedestalBegin;
        int pedEnd = pedestalEnd;
        int sigBegin = signalBegin;
        int sigEnd = signalEnd;


        for (int i = 0; i < 6; i++)
        {
          datafiles[iCh].read((char *)&buf, sizeof(buf));
        }
        /*
        if(iCh == 5 || iCh == 9){

          for (int i = 0; i < 6; i++)
          {
            datafiles[iCh].read((char *)&buf, sizeof(buf));
          }

        }
        */
        

        for (int i = 0; i < ev.recLength; i++)
        {
          datafiles[iCh].read((char *)&buf, sizeof(buf));
          int value = buf; // float to int here
          // if (value<0 || value>4095) value = 0;
          if (value < 0 || value > 8192)
            value = 0;
          ev.data[i] = (short)value;
        }

        //median filter
        std::vector<float> evdata(ev.data, ev.data + 1024);
        int filter_length = 5;
        for(int j=0; j<evdata.size()-filter_length; j++) {
          std::vector<short> sub(&evdata[j],&evdata[j+filter_length]);
          std::sort(sub.begin(), sub.end());
          ev.data[j] = sub.at(2);
        }      


        //Zero level calculation
        const int n_gates = 3;
        int gate_npoints = (int)floor((sigBegin-2.)/n_gates);

        float gates_mean[n_gates], gates_rms[n_gates];
        for(Int_t igate = 0; igate < n_gates; igate++)
          MeanRMScalc(ev.data, gates_mean+igate, gates_rms+igate, igate*gate_npoints, (igate+1)*gate_npoints);

        int best_gate = 0;
        for(Int_t igate = 0; igate < n_gates; igate++)
          if(gates_rms[igate] < gates_rms[best_gate]) best_gate = igate;

        ev.avgPed = gates_mean[best_gate];

        int time_max_in_gate_ = 0;
        for(UShort_t i = sigBegin; i <= sigEnd; i++)
        {
          float value = ev.data[i];
          if (signalType[iCh] < 0)
          { // falling
            if (value < ev.min){
              ev.min = value;
              time_max_in_gate_ = i;
            }
            dataInt += ev.avgPed - value;
            ev.dataInt[i] = dataInt;
          }
          else
          { // rising
            if (value > ev.max){
              ev.max = value;
              time_max_in_gate_ = i;
            }
            dataInt += value - ev.avgPed;
            ev.dataInt[i] = dataInt;
          }
        }
        ev.dataIntegral = dataInt;

        // calculate the time at threshold
        float threshValue = (ev.avgPed - ev.min) * constantFractionThresh; // constant fraction
        if (signalType[iCh] > 0)
        {
          threshValue = (ev.max - ev.avgPed) * constantFractionThresh;
        }
        //if(iCh == fit_channel) threshValue = 4096.0 / 2.0 * 0.01; // 10mV no cfd

        if (showData > 0)
          std::cout << "threshValue [" << iCh << "] = " << threshValue << std::endl;

        int point;
        point = time_max_in_gate_;
	      float front_time = (signalType[iCh] > 0) ? GoToLevel(ev.data, threshValue+ev.avgPed, &point, -1, sigEnd) :
                                                   GoToLevel(ev.data, ev.avgPed-threshValue, &point, -1, sigEnd) ;
        point = time_max_in_gate_;
	      float back_time = (signalType[iCh] > 0) ? GoToLevel(ev.data, threshValue+ev.avgPed, &point, 1, sigEnd) :
                                                  GoToLevel(ev.data, ev.avgPed-threshValue, &point, 1, sigEnd) ;

        
        if(iCh == com_channel1 || iCh == com_channel2) {
          short int comp_diff[ev.recLength];
          for (int i = 0; i < ev.recLength; i++)
          {
            if(iCh == com_channel1) datafiles[com_channel1+1].read((char *)&buf, sizeof(buf));
            if(iCh == com_channel2) datafiles[com_channel2+1].read((char *)&buf, sizeof(buf));
            int value = buf; // float to int here
            // if (value<0 || value>4095) value = 0;
            if (value < 0 || value > 8192)
              value = 0;
            comp_diff[i] = ev.data[i]-(short)value;
          }
          point = time_max_in_gate_;
          front_time = GoToLevel(comp_diff, 0, &point, -1, sigEnd);
          point = time_max_in_gate_;
	        back_time = GoToLevel(comp_diff, 0, &point, 1, sigEnd);
        }
        

        ev.timeOverThresh = back_time - front_time;
        ev.timeAtThresh = front_time;

        channelData.push_back(ev);

      } // for(int iCh = 0; iCh < N_CHANNELS; iCh++)

      // analyse data of one event here

      if (iEvt % 1000 == 0)
        cout << "Event: " << iEvt << endl;

      double adc[N_CHANNELS];
      double qdc[N_CHANNELS];
      double tdc[N_CHANNELS];
      double tot[N_CHANNELS];
      double ped[N_CHANNELS];

      for (int ch = 0; ch < N_CHANNELS; ch++)
      {
        adc[ch] = 0;
        qdc[ch] = 0;
        tdc[ch] = 0;
        tot[ch] = 0;
        ped[ch] = 0;
      }

      for (int ch = 0; ch < N_CHANNELS; ch++)
      {
        if (signalType[ch] < 0)
        { // falling
          adc[ch] = channelData.at(ch).avgPed - channelData.at(ch).min;
        }
        else
        {
          adc[ch] = channelData.at(ch).max - channelData.at(ch).avgPed;
        }
        qdc[ch] = channelData.at(ch).dataIntegral;
        tdc[ch] = channelData.at(ch).timeAtThresh;
        tot[ch] = channelData.at(ch).timeOverThresh;
        // ped[ch] = channelData.at(ch).max - channelData.at(ch).minPed;
        ped[ch] = channelData.at(ch).avgPed;

        _treeFill_adc[ch] = adc[ch];
        _treeFill_qdc[ch] = qdc[ch];
        _treeFill_tdc[ch] = tdc[ch];
        _treeFill_tot[ch] = tot[ch];
        _treeFill_ped[ch] = ped[ch];
      }

      _bigTree->Fill();

      //if(adc[fit_channel] > 3500) continue;
      std::vector<float> wfm(channelData.at(fit_channel).data, channelData.at(fit_channel).data + 1024);
      bool try_fit = !(wfm.empty() || tdc[fit_channel]<signalBegin || tdc[fit_channel]>signalEnd-20 || adc[fit_channel]>3500);
      float fit_integral, fit_chi2, fit_R2;
      std::vector<float> fit_wfm;
      if(try_fit) {

        PsdSignalFitting::PronyFitter Pfitter(model_order, exponents, signalBegin, signalEnd);
        //Pfitter.SetDebugMode(1);
        complex<float> *harmonics;
        Pfitter.SetWaveform(wfm, ped[fit_channel]);
        
        int best_signal_begin = Pfitter.ChooseBestSignalBeginHarmonics(tdc[fit_channel]-5, tdc[fit_channel]+20);
        Pfitter.SetSignalBegin(best_signal_begin);
        Pfitter.CalculateFitHarmonics();
        Pfitter.CalculateFitAmplitudes();
        
        /*
             std::vector<std::complex<float>> ext_harmonics = { {0.88,0}, {0.977,0}}; //laser 5ghz
            Pfitter.SetExternalHarmonics(ext_harmonics);
            int best_signal_begin = Pfitter.ChooseBestSignalBegin(tdc[fit_channel]-5, tdc[fit_channel]+20);
            Pfitter.SetSignalBegin(best_signal_begin);
            Pfitter.CalculateFitAmplitudes();
        */
        harmonics = Pfitter.GetHarmonics();
        auto amplituds = Pfitter.GetAmplitudes();
        fit_integral = Pfitter.GetIntegral(signalBegin, signalEnd);
        fit_chi2 = Pfitter.GetChiSquare(signalBegin, signalEnd, tdc[fit_channel]);
        fit_R2 = Pfitter.GetRSquare(signalBegin, signalEnd);
        if(false) printf("fit integral %.0f integral %.0f chi2 %.1f R2 %.3f\n", fit_integral, qdc[fit_channel], fit_chi2, fit_R2);
        fit_wfm = Pfitter.GetFitWfm();

        for(int wfm_iter=0; wfm_iter<1024; wfm_iter++){
          if(wfm_iter>=tdc[fit_channel]-10)
            wfm_ethalon->Fill(wfm_iter - tdc[fit_channel]+10, wfm.at(wfm_iter));
        }

        if(fit_R2<fit_R2_thr)
        {
          
            //Pfitter.SetDebugMode(1);
            //Pfitter.SetSignalBegin(best_signal_begin);
            //Pfitter.CalculateFitHarmonics();
            //Pfitter.CalculateFitAmplitudes();
          
            counter++;
            for(int i =0; i < exponents; i++)
                harmHistVect.at(i)->Fill(real(harmonics[i+1]), imag(harmonics[i+1]));

            for(int i =0; i < exponents; i++)
                tauHistVect.at(i)->Fill(-0.2/log(real(harmonics[i+1])));

            c_corr = 0;
            for(int i =0; i < exponents; i++)
              for(int j =0; j < exponents; j++){
                if(real(amplituds[1])>0) continue;
                amplCorrVect.at(c_corr)->Fill(0.488*real(amplituds[j+1]), 0.488*real(amplituds[i+1])); 
                c_corr++;
              }

            //cout<<"Harmonics"<<endl;
            //for(int i =0; i < exponents; i++)
            //    cout<<"H"<<i<<" "<< -0.2/log(real(harmonics[i+1]))<<" A "<<amplituds[i+1]<<"    ";
            //cout<<endl<<endl;

            //##################################
            TString signal_name = Form("channel %i fit_integral %.0f integral %.0f chi2 %.1f fit_R2 %.3f", 
                                fit_channel, fit_integral, qdc[fit_channel], fit_chi2, fit_R2);
            //##################################       
            //cout<<signal_name<<endl;             
                                
            std::vector<float> points(wfm.size());
            std::iota(std::begin(points), std::end(points), 0); // Fill with 0, 1, ..., wfm.back().
            if(!fit_wfm.empty()){
                TGraph *tgr_ptr_fit = new TGraph(fit_wfm.size(), &points[0], &fit_wfm[0]);
                tgr_ptr_fit->SetLineColor(kRed);
                tgr_ptr_fit->SetLineWidth(2);
                tgr_ptr_fit->Draw("same");
            }

        }
      }


      if (showData == 1 && (eventToShow < 0 || (eventToShow > 0 && iEvt == eventToShow)))
      {
        // if (isMuon) {

        if (eventToShow > 0 && iEvt == eventToShow)
          eventToShow = -1;

        cout << "iEvt = " << iEvt << endl;

        if (!(c1))
          c1 = new TCanvas("c1", "c1", 800, 600);
        c1->SetGrid();

        TGraph *gr[N_CHANNELS_TO_SHOW];
        for (int iGr = 0; iGr < N_CHANNELS_TO_SHOW; iGr++)
          gr[iGr] = new TGraph();

        for (int iCh = 0; iCh < N_CHANNELS_TO_SHOW; iCh++)
        {
          cout << "  ,ped [" << iCh << "] = " << ped[iCh];
        }
        cout << endl;

        for (int iCh = 0; iCh < N_CHANNELS_TO_SHOW; iCh++)
        {
          cout << "  ,adc [" << iCh << "] = " << adc[iCh];
        }
        cout << endl;

        for (int iCh = 0; iCh < N_CHANNELS_TO_SHOW; iCh++)
        {
          cout << "  ,qdc [" << iCh << "] = " << qdc[iCh];
        }
        cout << endl;

        for (int iCh = 0; iCh < N_CHANNELS_TO_SHOW; iCh++)
        {
          cout << "  ,tdc [" << iCh << "] = " << tdc[iCh];
          for (int iData = 0; iData < 1024; iData++)
          {
            gr[iCh]->SetPoint(iData, iData, channelData.at(iCh).data[iData]);
          }
        }

        cout << endl;

        for (int iCh = 0; iCh < N_CHANNELS_TO_SHOW; iCh++)
        {
          cout << "  ,tot [" << iCh << "] = " << tot[iCh];
        }
        cout << endl;

        if (!(hf1))
          hf1 = new TH1F("hf1", "hf1", 100, 0, 1024);

        hf1->SetStats(0);
        hf1->GetYaxis()->SetRangeUser(rangeUser1, rangeUser2);
        hf1->Draw();

        int colors[16] = {1, 2, 4, 6, 1, 2, 4, 6, 1, 1, 1, 1, 1, 1, 1, 1};

        for (int iGr = 0; iGr < N_CHANNELS_TO_SHOW; iGr++)
        {
          gr[iGr]->SetMarkerStyle(20);
          gr[iGr]->SetMarkerSize(0.5);
          gr[iGr]->SetMarkerColor(colors[iGr]);
          gr[iGr]->SetLineColor(colors[iGr]);
          gr[iGr]->Draw("plsame");
        }

        if(fit_R2<fit_R2_thr && try_fit)
        {
            //##################################
            TString signal_name = Form("channel %i fit_integral %.0f integral %.0f chi2 %.1f fit_R2 %.3f", 
                                fit_channel, fit_integral, qdc[fit_channel], fit_chi2, fit_R2);
            //##################################       
            cout<<signal_name<<endl;             
                                
            std::vector<float> points(wfm.size());
            std::iota(std::begin(points), std::end(points), 0); // Fill with 0, 1, ..., wfm.back().
            if(!fit_wfm.empty()){
                TGraph *tgr_ptr_fit = new TGraph(fit_wfm.size(), &points[0], &fit_wfm[0]);
                tgr_ptr_fit->SetLineColor(kRed);
                tgr_ptr_fit->SetLineWidth(2);
                tgr_ptr_fit->Draw("same");
            }
        }



        TLine *line1 = new TLine(0, ped[0], 1024, ped[0]);
        line1->SetLineColor(2);
        line1->Draw();

        Double_t rangeUser3 = rangeUser1 + (rangeUser2 - rangeUser1) * 0.2;
        Double_t rangeUser4 = rangeUser2 - (rangeUser2 - rangeUser1) * 0.2;

        TLine *line2 = new TLine(pedestalBegin, rangeUser3, pedestalBegin, rangeUser4);
        TLine *line3 = new TLine(pedestalEnd, rangeUser3, pedestalEnd, rangeUser4);
        line2->SetLineColor(2);
        line2->SetLineWidth(2);
        line3->SetLineColor(2);
        line3->SetLineWidth(2);
        line2->Draw();
        line3->Draw();

        TLine *line4 = new TLine(signalBegin, rangeUser3, signalBegin, rangeUser4);
        TLine *line5 = new TLine(signalEnd, rangeUser3, signalEnd, rangeUser4);
        line4->SetLineColor(2);
        line4->SetLineWidth(2);
        line5->SetLineColor(2);
        line5->SetLineWidth(2);
        line4->Draw();
        line5->Draw();

        c1->Update();
        char c = getchar();
        if (c == 'q')
          break;
        else if (c == 'c')
          showData = 0;
        else if (c == 'p')
          c1->SaveAs(Form("c1_%d.png", iEvt));
      }

    } // for (int iEvt=0; iEvt<nEvents; iEvt++)

    _bigTree->Write();

    outfile.Close();

    TCanvas *canv_wfm = new TCanvas();
    canvas_array.push_back(canv_wfm);
    wfm_ethalon->Draw();

    for(int i =0; i < exponents; i++){
        TCanvas *canv_ptr = new TCanvas();
        canvas_array.push_back(canv_ptr);
        tauHistVect.at(i)->Draw("colz");
    }

    TCanvas *canv_taus = new TCanvas();
    canvas_array.push_back(canv_taus);
    for(int i =0; i < exponents; i++){
        tauHistVect.at(i)->SetLineColor(40+i);
        tauHistVect.at(i)->Draw("same");
    }

    TCanvas *canv_corr = new TCanvas();
    canvas_array.push_back(canv_corr);
    amplCorrVect.at(1)->Draw("colz"); 
    /*
    canv_corr->Divide(exponents, exponents);
    c_corr = 0;
    for(int i =0; i < exponents; i++)
      for(int j =0; j < exponents; j++){
        canv_corr->cd(c_corr+1);
        amplCorrVect.at(c_corr)->Draw("colz"); 
        c_corr++;
      }
    */

    for (Int_t i = 0; i < canvas_array.size() - 1; i++)
      ((TCanvas *)canvas_array.at(i))->SaveAs((pics_file_name + ".pdf(").Data());
    ((TCanvas *)canvas_array.back())->SaveAs((pics_file_name + ".pdf)").Data());
    printf("file %s was written\n", (pics_file_name + ".pdf").Data());

    saveHistogramsToRootFile(canvas_array, (pics_file_name + ".root").Data());

    cout<<"total good fits  "<<counter<<endl;

  }
  else
  {
    cout << "Unable to open input files" << endl;
  }
}
