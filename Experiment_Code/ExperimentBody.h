#pragma once

#include "emil/EMIL.hpp"
#include <fstream>
#include "Pest.h"

class BadConversion : 
public std::runtime_error 
{
public:
BadConversion(std::string const& s)
: std::runtime_error(s)
{ }
};
inline std::string stringify(double x)
 {
  std::ostringstream o;
 if (!(o << x))
  throw BadConversion("stringify(double)");
return o.str();
} 

class ExperimentBody : public CExperiment
{
public:
	ExperimentBody(int pxWidth, int pxHeight,
		int RefreshRate, CCfgFile* Params);

	/// Standard event handlers
	void initialize();
	void finalize();
	void eventRender(unsigned int FrameCount, CEOSData* Samples);

	void eventJoypad();
	void eventKeyboard(unsigned char key, int x, int y);

	void createTrialSequence();
	float getGazePxDist(float pX, float pY); // get distance between gaze to input point
	//////////////////////

private:

	enum STATE
	{
		STATE_LOADING,
		STATE_CALIBRATION,
		STATE_PERIPH_FIX,
		STATE_SUSTAIN_FIX,
		STATE_CENTRAL_FIX,
		STATE_IMAGE,
		STATE_WAIT_RESPONSE,
		STATE_QUIT
	};

	STATE m_state;


	void setParameters();
	void createObjects();
	void saveData();

	// Configuration file
	CCfgFile* m_paramsFile;

	string m_sbj;

	// Stimuli and parameters for the test calibration
	CSolidPlane* m_fixcross;
	CSolidPlane* m_gazecross;
	int Increment;

	// fixation stimuli - central and peripheral cues
	CSolidPlane* m_pCue;
	CSolidPlane* m_cCue;
	CImagePlane* m_crosshair;

	// shaders and image objects for gratings
	std::shared_ptr<CShader> m_imageShader;
	CImagePlane* m_noiseBackground;

	// timers
	CTimer m_trialTimer; // full trial
	CTimer m_stimTimer; // stimulus presentation
	CTimer m_rt; // response time

	float ResponseTime;

	// experiment parameters
	int m_trialNo;
	int m_rampTime;
	int respWaitTime;
	std::string gratingType;
	int m_spatialFreqLev;
	int m_eccLevUsed;
	float pixelAngle;
	int m_backLum;
	float m_contrastBack;
	//float eccentricity;

	int m_currTrial;

	// temporary trial variables
	struct trialParams {
		float contrast; // varies with pest if requested
		// variables that definitely vary trial to trial
		float rt;
		float fixTime;
		int resp;
		int seed;//to store value of random_device
		int presTime;
		int presLev;
		int ecc;
		int spFreq;
		int instance;		
		float fixAngle;
		float phase;
		int present;
		string background_img;
		string central_img;
	};

	vector <trialParams> m_trial;
	trialParams curr_trial;

	// screen locations
	float centerX;
	float centerY;
	float periphX;
	float periphY;
	float landingX;
	float landingY;

	// time markers
	float cuePrf;
	float fixOn;
	float cueCtr;
	float saccOn;
	float saccOff;
	float flashOn;
	float rampOff;
	float stimOff;
	float quit;

	float m_pestStep;
	float m_pestTarget;
	float m_saccDur;
	float m_fixEcc;
	float m_fixAngleOffset;
	float m_fixAreaC;
	float m_fixAreaP;
	float m_fixAreaCpx;
	float m_fixAreaPpx;
	

	// general stimulus information
	int rampFrame;
	int currRampFrame;
	int currFrame;
	float awidth;
	float gaussSD;
	float gaussSD_2;
	
	int m_numCompleted;
	int m_numTestCalibration;
	int m_numSession;
	bool debug;

	int Response;

	float m_presTime[2];
	float m_contrast[7];
	float m_ecc[3];
	float m_spFreq[2];


	float m_respTime;
	float m_fixTime;
	float m_fixTimeRng;
	float m_loadTime;
	float m_dropOutTime;
	bool m_stabilized;
	bool m_flashed;
	bool m_fixMode;
	bool m_feedback;

	int m_presTimeLevUsed;
	int m_contrastLevUsed;
	int m_totalInstances;
	int m_fixLocNo;
	int m_imgNo;

	// eye position information
	float Dist;
	float X; // raw eye position (px)
	float Y;
	float Xstab; // stab eye position (px)
	float Ystab;
	float xpos; // stimulus position
	float ypos;
	float m_corrX;
	float m_corrY;

	bool m_buttonPress;
	bool m_isOpen;

	bool SHOW_TARGET;
	float gain;
	int nRecal;
	float refreshRate;

	vector <Pest*> m_pest;
	//vector <vector<Pest*>> m_pest;
};