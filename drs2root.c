//
//  drs2root.C
//
//  convert DRS4 data to root files
//

#include <iostream>
#include <stdio.h>
// #include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <cmath>
#include <TH2.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>

#include "Compression.h" 

using namespace std;

const int NSAMPLES=1024;

typedef struct {
   char time_header[4]; // 'T' 'I' 'M' 'E'
   char bn[2];          // 'B' '#'
   unsigned short board_number;
} THEADER;

typedef struct {
   char channel[4];     // eg, 'C' '0' '0' '1'
   float tcal[1024];
} TCHEADER;

typedef struct {
   char event_header[4];
   unsigned int ev_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisec;
   unsigned short range;
   char bn[2];
   unsigned short board_number;
   char tc[2];
   unsigned short trigger_cell;
} EHEADER;

typedef struct {
   char channel[4];
   unsigned short data[NSAMPLES];
} CHANNEL;

typedef struct {
   char channel[4];
   unsigned int scaler;  // new?, was present in output from drs-5.0.5)
   unsigned short data[NSAMPLES];
} CHANNEL_WSCALER;


/*-----------------------------------------------------------------------------*/

#define N_DEVCHN 4 // number of channels in device
#define N_CHN 4 // number of saved channels


int drs2root(const char *fname = "Finger_azb_nov_4_2018_70cm.dat", const int do_display = 0)
{
   int ch, i, j, fh, ndt;
   double threshold, t1, t2, dt, sumdt, sumdt2;
   THEADER th;
   TCHEADER tch[N_CHN];
   EHEADER eh;
   CHANNEL_WSCALER channel[N_CHN];
   //double waveform[N_CHN][NSAMPLES], time[N_CHN][NSAMPLES];
   float waveform[N_CHN][NSAMPLES], time[N_CHN][NSAMPLES];
   
   const int do_check = 0;
   int nevents = 0;
   int n4hitevents = 0;
   int nhits[N_CHN] = {0};  // how many above threshold for each ch

   // open file
   fh = open(fname, O_RDONLY, 0666);
   if (fh < 0) {
      printf("Cannot find data file\n");
      return 0;
   }

   // read time header (one set of calibs per file)
   read(fh, &th, sizeof(th));
   if ( do_check ) {
     cout << "board no. " << th.board_number << endl;
   }
   // read time bin widths (tcal)
   for (ch=0 ; ch<N_CHN ; ch++) {
     read(fh, &tch[ch], sizeof(tch[ch]));
     if ( do_check ) {
       cout << "====  channel " << ch << endl;
       for (int i=0; i<NSAMPLES; i++)
       {
         cout << tch[ch].tcal[i] << "\t";
         if ( i%8==7 ) cout << endl;
         tch[ch].tcal[i] = 1.5;
       }
     }
   }

   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;

   TCanvas *c_display = new TCanvas("c_display","FTBF 2016",1000,400);
   TGraph *g_trace[4];
   TString name;
   for (int ich=0; ich<4; ich++) {
     name = "ch"; name += (ich+1);
     g_trace[ich] = new TGraph();
     g_trace[ich]->SetName(name);
   }

   // open the output root file
   TString savefname = fname;
   savefname.ReplaceAll(".dat",".root");
   //TFile *outfile = new TFile(savefname, "RECREATE", "", 9);
   TFile *outfile = new TFile(savefname, "RECREATE", "", ROOT::kLZMA*100 + 9);
   //TFile *outfile = new TFile(savefname, "RECREATE", "", 0);
   
   // define the rec tree
   Int_t f_evt;
   Int_t f_spill;
   Int_t f_spillevt;
   TTree *rec = new TTree("t","DRS4");
   rec->Branch("evt", &f_evt,"evt/I");  
   rec->Branch("spill", &f_spill,"spill/S");  
   rec->Branch("spillevt", &f_spillevt,"spillevt/S");  
/*
   rec->Branch("t0", &time[0]   ,"t0[1024]/D");  
   rec->Branch("t1", &time[1]   ,"t1[1024]/D");  
   rec->Branch("t2", &time[2]   ,"t2[1024]/D");  
   rec->Branch("t3", &time[3]   ,"t3[1024]/D");  
   rec->Branch("ch0", &waveform[0] ,"ch0[1024]/D");
   rec->Branch("ch1", &waveform[1] ,"ch1[1024]/D");
   rec->Branch("ch2", &waveform[2] ,"ch2[1024]/D");
   rec->Branch("ch3", &waveform[3] ,"ch3[1024]/D");
   */
   // changed to floats to save some memory and disk space
   rec->Branch("t0", &time[0]   ,"t0[1024]/F");  
   rec->Branch("t1", &time[1]   ,"t1[1024]/F");  
   rec->Branch("t2", &time[2]   ,"t2[1024]/F");  
   rec->Branch("t3", &time[3]   ,"t3[1024]/F");  
   rec->Branch("ch0", &waveform[0] ,"ch0[1024]/F");
   rec->Branch("ch1", &waveform[1] ,"ch1[1024]/F");
   rec->Branch("ch2", &waveform[2] ,"ch2[1024]/F");
   rec->Branch("ch3", &waveform[3] ,"ch3[1024]/F");


  
   /*
   rec = new TTree("t","mcp-pmt data");
   rec->Branch("evt",&f_evt,"evt/I");
   rec->Branch("tstamp",&f_tstamp,"tstamp/I");
   rec->Branch("spill",&f_spill,"spill/S");
   rec->Branch("spillevt",&f_spillevt,"spillevt/S");
   rec->Branch("ch",&f_ch,"ch[12]/S");
   rec->Branch("v",f_volt,"v[12][1024]/F");    // waveform value in volts
   rec->Branch("t",f_time,"t[12][1024]/F");    // waveform time in ns
   rec->Branch("ampl",&f_ampl,"ampl[12]/F");
   rec->Branch("q",&f_q,"q[12]/F");
   rec->Branch("t0",&f_t0,"t0[12]/F");
   */

   

   int hit[N_CHN] = {0};  // whether all 4 were hit
   int nhit = 0;  // how many above threshold
   for (int nevt= 0 ;nevt<1000000 ; nevt++) {
     
      bool looking_for_next_evt = true;
      while (looking_for_next_evt) {
        memset(&eh, sizeof(eh), 0);
        // read event header
        i = (int)read(fh, &eh, sizeof(eh));

        // check for valid event header
        if (memcmp(eh.event_header, "EHDR", 4) != 0) {
          //printf("Invalid event header (probably number of saved channels not equal %d)\n", N_CHN);
          printf("Invalid event header #%d.  Searching for next event header...\n",nevt);
          lseek(fh, off_t(-sizeof(eh) + 4), SEEK_CUR); // try to step forward 1 word (4 bytes) before trying to read the event header
        } else {
          looking_for_next_evt = false;
        }

      }
      f_evt = nevt+1;

      if (do_check)
      {
        cout << "evt\t" << eh.ev_serial_number << endl;
        cout << "range\t" << eh.range << endl;
        cout << "board\t" << eh.board_number << endl;
        cout << "trigcell\t" << eh.trigger_cell << endl;
 
        cout << "Date/Time = " << eh.year << eh.month << eh.day <<  " " << eh.hour << ":"<< eh.minute << ":"<< eh.second << "." << eh.millisec << endl;

      }
// return 0;
      //cout << "Range = " << eh.range << endl;

      // print notification every 100 events
      if (nevt % 100 == 0 || i != sizeof(eh))
         printf("Analyzing event #%d\n", nevt);
      
      // stop if end-of-file
      if (i != sizeof(eh))
         break;
      
      nhit = 0;
      // read channel data
      for (int ch=0 ; ch<N_CHN ; ch++)
      {

        //cout << "ch " << ch << endl;
        hit[ch] = 0;

        read(fh, &channel[ch], sizeof(channel[ch]));

        // check for valid channel header
        if (memcmp(channel[ch].channel, "C00", 3) != 0) {
          printf("Invalid channel header  %d\n", ch);
          return 0;
        }

        for (int i=0 ; i<NSAMPLES ; i++)
        {
          // convert to volts
          // this is the older version
          //waveform[ch][i] = (channel[ch].data[i]/65536.0-0.5)*1000.;
          // this is what is used in the latest code (v5.0.4)
          waveform[ch][i] = (channel[ch].data[i] / 65536. + eh.range/1000.0 - 0.5);

          // calculate time for this cell
          time[ch][i] = 0;
          for (int j=0; j<i ; j++)
          {
            time[ch][i] += tch[ch].tcal[(j+eh.trigger_cell) % NSAMPLES];
          }
        }
      }

      //cout << "Aligning cells" << endl;

      // align cell #0 of all channels
      t1 = time[0][(NSAMPLES-eh.trigger_cell) % NSAMPLES];
      for (int ch=1 ; ch<N_CHN ; ch++) {
        t2 = time[ch][(NSAMPLES-eh.trigger_cell) % NSAMPLES];
        dt = t1 - t2;
        for (int i=0 ; i<NSAMPLES ; i++)
        {
          time[ch][i] += dt;
        }
      }

      //cout << "writing tree" << endl;

      rec->Fill();

      //cout << "done writing tree" << endl;
      // Find how many are above threshold
      for (int ch=0 ; ch<N_CHN ; ch++)
      {
        for (int i=0 ; i<NSAMPLES ; i++)
        {
          if (time[ch][i]>10.&&time[ch][i]<40.&&waveform[ch][i]>0.02)
          {
            hit[ch] = 1;
            break;
          }
        }

        if ( hit[ch]==1 )
        {
          nhit++;
          nhits[ch]++;
          for (int i=NSAMPLES/2 ; i<NSAMPLES ; i++)
          {
            // veto on noise events from current fluctuations
            if (time[ch][i]>130.&&time[ch][i]<150.&&waveform[ch][i]>0.1)
            {
              nhit--;
              nhits[ch]--;
              break;
            }
          }
        }
      }

      // draw every do_display events
      if ( do_display!=0 && nevt%do_display==0 && hit[0]==1 ) {
      //if ( do_display!=0 && nhit==4 ) {
        rec->SetLineColor(1);
        rec->Draw("ch0:t0","evt==1","pl");
        rec->SetLineColor(2);
        rec->Draw("ch1:t1","evt==1","plsame");
        rec->SetLineColor(3);
        rec->Draw("ch2:t2","evt==1","plsame");
        rec->SetLineColor(4);
        rec->Draw("ch3:t3","evt==1","plsame");

        //gPad->ls();
        TH2 *htemp = (TH2*)gPad->GetPrimitive("htemp");
        if ( htemp!=0 )
        {
          htemp->GetYaxis()->SetRangeUser(-0.05,0.2);
          htemp->GetXaxis()->SetRangeUser(10.,40.);
        }
        gPad->Modified();
        gPad->Update();
        string junk;
        cin >> junk;
        if (junk[0] == 'q') break;
      }

      nevents++;
      if ( nhit==4 )
      {
        //cout << "found hit=4" << endl;
        n4hitevents++;
      }

   }
 
   // save and close root file
   rec->Write();
   outfile->Write();
   //outfile->Close();

   cout << "Processed " << nevents << endl;
   if ( nevents!=0 && n4hitevents!=0 )
   {
     double fraction = (double)n4hitevents/nevents;
     double dfraction = fraction*sqrt(1.0/n4hitevents + 1.0/nevents);
     cout << "Found " << n4hitevents << " with hits=4\t" << fraction << " +- " << dfraction << endl;
   }
   else
   {
     cout << "Found " << n4hitevents << " with hits=4\t" << endl;
   }

   // for back mrpc
   if ( nevents!=0 && nhits[2]!=0 )
   {
     double fraction = (double)nhits[2]/nevents;
     double dfraction = fraction*sqrt(1.0/nhits[2] + 1.0/nevents);
     cout << "Found " << nhits[2] << " with hit in back mRPC\t" << fraction << " +- " << dfraction << endl;
   }

   cout << "hits per channel" << endl;
   for (int ich=0; ich<N_CHN; ich++)
   {
     cout << ich << "\t" << nhits[ich] << endl;
   }

   return 1;
}

