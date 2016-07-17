#include <cstdlib>
#include <iostream>
#include <cstring>
#include <time.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <cstring>
#include <string>
#include "ErrorFinderManager.hpp"
#include "ErrorCalculator.hpp"

#include "Consolidator.hpp"



#define GetCurrentDir getcwd

/*This function was copied from consolidator.cpp. This is a standalone function*/
template <typename T>
  std::string NumberToString ( T Number )
  {
     std::ostringstream ss;
     ss << Number;
     return ss.str();
  }

ErrorFinderManager::ErrorFinderManager():WINDOW(50),
                    MIN_SNP(30),GAP(1), MA_SNP_END(50), TRUESNP( 600 ),
                    MIN_CM(0.4), MA_ERR_THRESHOLD_START(0.08),
                    MA_ERR_THRESHOLD_END(0.08),PCT_ERR_THRESHOLD( 0.90 ),
                    HO_THRESHOLD( 0.98 ), TRUECM( 6 ),PIELENGTH( 3 ),
                    ISMOL( false), COUNTGAPERR( false ),MA_THRESHOLD(0.8),EMPIRICAL_MA_RESULT(-1.0),EMPIRICAL_PIE_RESULT(-1.0) 
{
}

void ErrorFinderManager::performConsolidation(int argc,char *argv[])
{
       bool goodParam=true;
       bool thresholdError = false; //if the user supplies both -empricial-ma-threshold and -ma-threshold, that is an error. This will be used to detect such an error.
       bool pieThresholdError = false;
       char currentPath[FILENAME_MAX];
       if(!GetCurrentDir(currentPath,sizeof(currentPath))){
std::cerr << "Error reading current directory" << std::endl;
		return;
       }
       for(int i=1;i<argc;i++)
        {

                if(strcmp(argv[i],"-bmatch")==0&&i<argc-1)
                {
                     BMATCHFILE=std::string(argv[++i]);
                }
                else if(strcmp(argv[i],"-bmid")==0&&i<argc-1)
                {
                     BMIDFILE=std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-bsid")==0&&i<argc-1)
                {
                     BSIDFILE=std::string(argv[++i]);
                }
                 else if(strcmp(argv[i],"-reduced")==0&&i<argc-2)
                {
                     MIN_SNP=atoi(argv[++i]);
                     MIN_CM=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ped-file")==0&&i<argc-1)
                {
                     PEDFILE=std::string(argv[++i]);
                }
		            else  if(strcmp(argv[i],"-holdout-ped")==0&&i<argc-1)
                {
                     HPEDFILE=std::string(argv[++i]);
                }
		            else  if(strcmp(argv[i],"-holdout-map")==0&&i<argc-1)
                {
                     HMAPFILE=std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-window")==0&&i<argc-1)
                {
                     WINDOW=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-err-threshold-start")==0&&i<argc-1)
                {
                     MA_ERR_THRESHOLD_START=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-holdout-threshold")==0&&i<argc-1)
                {
                     HO_THRESHOLD=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-trueCM")==0&&i<argc-1)
                {
                     TRUECM=atof(argv[++i]);
                }
                else  if( strcmp( argv[i],"-trueSNP" )==0 && i < argc-1 )
                {
                     TRUESNP=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-holdout-missing")==0&&i<argc-1)
                {
                     HO_MISSING= std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-err-threshold-end")==0&&i<argc-1)
                {
                     MA_ERR_THRESHOLD_END=atof(argv[++i]);
                }
                else  if(strcmp(argv[i],"-gap")==0&&i<argc-1)
                {
                     GAP=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-ma-snp")==0&&i<argc-1)
                {
                     MA_SNP_END=atoi(argv[++i]);
                }
                else  if(strcmp(argv[i],"-pct-err-threshold")==0&&i<argc-1)
                {
                     if(pieThresholdError == true){
std::cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << std::endl;
                      exit(1);
                     }
                     PCT_ERR_THRESHOLD=atof(argv[++i]);
                     pieThresholdError = true;
                }
                else if(strcmp(argv[i],"-emp-pie-threshold")==0&&i<argc-1)
                {
                     if(pieThresholdError == true){
                    	 std::cerr << "ERROR: You have supplied both -emp-pie-threshold and -pct-err-threshold parameters, but only one is allowed. Please try again." << std::endl;
                      exit(1);
                     }
                     EMPIRICAL_PIE_RESULT=atof(argv[++i]); 
                     pieThresholdError = true;               
                }
		            else  if(strcmp(argv[i],"-output.type")==0&&i<argc-1)
                {
                     OPTION=std::string(argv[++i]);
                }
                else  if(strcmp(argv[i],"-log.file")==0&&i<argc-1)
                {
                     LOGFILE=std::string(argv[++i]);
                } 
                else if(strcmp(argv[i], "-trueibd")==0&&i<argc-1)
                {
                     TRUEIBDFILE = std::string(argv[++i]);
                }
            		else if(strcmp(argv[i],"-ma-threshold")==0&&i<argc-1)//adding new -ma-threshold argument
            		{
            		     if(thresholdError == true){ //user has already supplied an empirical-ma-threshold, so exit the program with an error message
            		    	 std::cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< std::endl;
                       exit(1);
            		     }
            		     MA_THRESHOLD=atof(argv[++i]);
            		     thresholdError = true;
            		}
            		else if(strcmp(argv[i],"-empirical-ma-threshold")==0 && i<argc-1){
            		     if(thresholdError == true){
            		    	 std::cerr << "ERROR: You have supplied both -empirical-ma-threshold and -ma-threshold parameters, but only one is allowed. Exiting program."<< std::endl;
            			exit(1);
            		     }	
            		     EMPIRICAL_MA_RESULT = atof(argv[++i]); //use the user supplied empirical ma threshold, instead of calculating it via true ibd segments	
            		     thresholdError = true;
            		}
                else  if(strcmp(argv[i],"-PIE.dist.length")==0&&i<argc-1)
                {
                     std::string MOL=std::string(argv[++i]);
                     if( MOL.compare( "MOL" ) ==0 )
                     {
                        ISMOL = true;
                     }
                     else
                     {
                        PIELENGTH = atof( MOL.c_str() );
                     }
                }
                else  if(strcmp(argv[i],"-count.gap.errors")==0&&i<argc-1)
                {
                     std::string option=std::string(argv[++i]);
                     if( option.compare( "TRUE" ) ==0 )
                     {
                        COUNTGAPERR = true;
                     }
                }

                else
                {
                        wrongParam += " " + std::string(argv[i]);

                        goodParam=false;
                }
        }
        if((!goodParam)||BMATCHFILE.compare("")==0||BSIDFILE.compare("")==0||BMIDFILE.compare("")==0||PEDFILE.compare("")==0)
        {

                displayError( argv[0] );
                return;
             
        }
        if( OPTION.compare( "" ) == 0 )
        {
        	std::cerr<< " please provide a valid output.type option " <<std::endl;
            exit( -1 );
        }
        if( LOGFILE.compare( "" ) == 0 )
        {
        	std::cerr<< " default log file name is FISH " <<std::endl;
            LOGFILE = "FISH";
        }
        eCalculator.createLogFile( LOGFILE  );
        eCalculator.countGapErrors( COUNTGAPERR ); 
        time_t startTime;

        time (&startTime);
        std::string str1 = " The program started at: " + std::string( ctime ( &startTime ) );
        std::string str = " Program working directory was: " + std::string(currentPath) +
		     " \nProgram version was: " + std::string(argv[0]) +
		     " \nProgram options:\n-bmatch file: " + BMATCHFILE +
                     " \n-bmid file: " + BMIDFILE +
                     " \n-bsid file: " + BSIDFILE +
                     " \n-ped file: " + PEDFILE +
                     " \n-holdout ped file: " + HPEDFILE +
                     " \n-holdoutmap file: " + HMAPFILE +
                     " \n-output type: " + OPTION +
                     " \n- missing SNP representation in pedfile: " + HO_MISSING +
                     " \n-log file: " + LOGFILE;
        str = str + " \nmin snp length : " + NumberToString( MIN_SNP )  +
                    " \nmin cm length : " + NumberToString( MIN_CM  ) +
                    " \ngap to consolidate : " + NumberToString( GAP ) +
                    " \nmoving averages window size : " + NumberToString(  WINDOW ) +
                    " \ndiscard ends to calculate pct err : " + NumberToString( MA_SNP_END ) +
                    "  \npercentage error threshold: " + NumberToString( PCT_ERR_THRESHOLD );
                  eCalculator.log( str );  

      initiateErrorFinder();
      //May need to stream in bmatch results at this point, or within the performTrim calculation itself.
      consolidator.setPersonCount(eCalculator.getNoOfPersons());
      if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
      {
    	  std::cerr<<"entering into HoldOut Mode Matching"<<std::endl;
           //updated version 5/30
	         consolidator.performTrim(eCalculator, WINDOW, MA_SNP_END,MA_THRESHOLD,MIN_SNP,MIN_CM,PCT_ERR_THRESHOLD,OPTION,HO_THRESHOLD,true,EMPIRICAL_MA_RESULT, BMATCHFILE, EMPIRICAL_PIE_RESULT);
	         std::cerr<<" Main Trim operation has completed "<< std::endl;
	         std::cerr<< " Hold out trim has completed" <<std::endl;
           if( (OPTION.compare( "finalOutput" ) == 0) || (OPTION.compare( "Full" ) == 0 ) )
           {   
                consolidator.finalOutPut( eCalculator, MIN_CM, MIN_SNP );
           }

      }
      else 
      {
          if( OPTION.compare( "Error3") == 0 )
          {
        	  std::cerr<< " Error: You have provided option:Error3 " << std::endl
                   <<" you can use Error3 only if you provided"
                   <<" hold out ped and map file, program with not output anything" << std::endl;
              exit( -1 );
          }
         
          consolidator.performTrim(eCalculator, WINDOW, MA_SNP_END, MA_THRESHOLD, MIN_SNP, MIN_CM, PCT_ERR_THRESHOLD, OPTION, HO_THRESHOLD, false,EMPIRICAL_MA_RESULT, BMATCHFILE, EMPIRICAL_PIE_RESULT);
		      
      }
      
      time_t endTime;
     
      time (&endTime);
      
      std::string str2 = " The program ended at: " + std::string( ctime ( &endTime ) );
      str2 = str2 +  "  Total time ( in seconds): " + NumberToString( ( endTime - startTime ) );
      eCalculator.log( str1 );
      eCalculator.log( str2 );    
}
void ErrorFinderManager::displayError(std::string argv)
{
	std::cerr<<"these parameters are not allowed "<<wrongParam<<std::endl;
	std::cerr << "Usage: " << argv << " -bmatch [BMATCH FILE]  -bsid [BSID FILE] -bmid [BMID FILE] -reduced [min_snp] [min_cm] "
  <<" -ped-file [ped file] -window [window width to calculate moving averages] "
  <<" -gap [max gap to consolidate two matches]"
  <<" -pct-err-threshold [max percentage of errors in a match after the trim] OR -emp-pie-threshold" 
  <<" -ma-threshold [specifies percentile to be drawn from trulyIBD data for MA calculations] OR -empirical-ma-threshold"
  <<" Note that if both -emp-pie-threshold and empirical-ma-threshold are supplied, then -trueSNP and -trueCM will be ignored"
  <<"-output.type [ must provide any of these. it can be "
  << "MovingAverages  or Error1 or Error2 or Error3 or ErrorRandom1 " 
  << "or ErrorRandom2 or Error3 or ErrorRandom3 or Full "
  <<  "look at the description about how these works in wiki ]"
  << "(optional) -holdout-ped [new ped file path] -holdout-map [new map file] "
  << "-holdout-threshold [threshold to drop a match with new ped file ]"
  << " -holdout-missing [missing value representation in new ped file] "
  << " -log.file [log file name]"
  << " -trueCM [ true match maximum cm length] " 
  << " - trueSNP [ true match SNP length]"
  << " -PIE.dist.length [ can be MOL or any cm distance length "
  << "please refer wiki for more details on how to use this option"
  << "-count.gap.errors [ TRUE or FALSE to include gap errors in errors count ]"
  << std::endl;
}
void ErrorFinderManager::initiateErrorFinder()
{
        eCalculator.readBmidFile(BMIDFILE);
        std::cerr<<"Reading bmid file completed"<<std::endl;
        eCalculator.readBsidFile(BSIDFILE);
        std::cerr<<"Reading bsid file completed"<<std::endl;
        eCalculator.readPedFile(PEDFILE, HO_MISSING);
        std::cerr<<"Reading ped file completed"<<std::endl;
        int pers_count=eCalculator.getNoOfPersons();
        if( HPEDFILE.compare( "" ) !=0 && HMAPFILE.compare( "" ) != 0 )
        {
          eCalculator.changeMapFile( HMAPFILE );
          eCalculator.readHPedFile( HPEDFILE, HO_MISSING );
          std::cerr<< " new map and ped File has been read" <<std::endl;
          std::cerr<< " calculating true percentage errors" <<std::endl;
          if( ISMOL )
          {
            consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count, TRUESNP, TRUECM);
          //    consolidator.findTruePctErrors( eCalculator, MA_SNP_END, true, WINDOW,MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          else
          {
            consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
            //  consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, true, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          std::cerr<< " true percentage errors calculated "<<std::endl;
       }
       else
       {
    	   std::cerr<< " calculating true percentage errors" <<std::endl;
          if( ISMOL )
          {
            consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
            //consolidator.findTruePctErrors( eCalculator, MA_SNP_END, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          else
          {
            consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT, BMATCHFILE, pers_count,TRUESNP, TRUECM );
            //consolidator.findTrueSimplePctErrors( eCalculator, PIELENGTH, false, WINDOW, MA_THRESHOLD, EMPIRICAL_MA_RESULT);
          }
          std::cerr<< " true hold out percentage errors calculated "<<std::endl;
      }   
  
}
