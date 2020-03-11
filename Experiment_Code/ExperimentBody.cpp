#include "stdafx.h"
#include "ExperimentBody.h"
#include "PostSaccSens_Ecc.h"
#include <math.h>
#include <iostream>
#include <random>

#define PI 3.14159265

/*
To Do:

To Check:

*/

///////////////////////////////////////////////////////////////////////////////////
ExperimentBody::ExperimentBody(int pxWidth, int pxHeight, int RefreshRate, CCfgFile* Params) :
	CExperiment(pxWidth, pxHeight, RefreshRate),

	m_paramsFile(Params)

{
	setExperimentName("PostSaccSens_Ecc");
	m_paramsFile = Params;
	CMath::srand(time(NULL));

	refreshRate = RefreshRate;
}

///////////////////////////////////////////////////////////////////////////////////
void ExperimentBody::initialize()
{

	CStabilizer::Instance()->enableSlowStabilization(false);


	COGLEngine::Instance()->clearScreen();
	glColor3d(255, 255, 255);
	//printCentered(CFontEngine::FONTS_ARIAL_18, 0, 0, "Loading...");

	// set parameters for online EM detection - this will likely change for each subject
	CEOS::Instance()->EMTM_setParameters(EOS_EYE_1, 15, 1000, 500, 5, 150, true, 30, 10);


	setParameters();
	createObjects();
	createTrialSequence();

	disable(CExperiment::EIS_PHOTOCELL);
	disable(CExperiment::EIS_NOTRACK_ICON);
	enable(CExperiment::EIS_STAT1);
	disable(CExperiment::EIS_JP_STRT);

	hideAllObjects();

	if (m_stabilized)
		CStabilizer::Instance()->enableSlowStabilization(true);

	// Set the background color for the experiment
	COGLEngine::Instance()->setBackgroundColor(m_backLum, m_backLum, m_backLum);

	/* Seed the random-number generator with current time so that
	* the numbers will be different every time we run.*/
	srand((unsigned)time(NULL));

	m_state = STATE_LOADING;
}




///////////////////////////////////////////////////////////////////////////////////
void ExperimentBody::eventRender(unsigned int FrameCount, CEOSData* Samples)
{

	// store the stabilized trace ( --- all in px)
	if (~debug)
	{
		X = m_corrX + Samples->x1;
		Y = m_corrY + Samples->y1;
	}
	else
	{
		X = 0 + m_corrX; Y = 0 + m_corrY;
	}

	// If the trial is running too long, just quit
	if (!m_fixMode && m_trialTimer.isExpired() &&
		m_trialTimer.isRunning()) {

		m_trialTimer.stop();
		m_state = STATE_QUIT;
	}

	float x;
	float y;

	switch (m_state) {
	case STATE_CALIBRATION:
		if (m_isOpen) { // start the calibration trial
			startTrial();
			m_gazecross->show();
			m_fixcross->show();
			moveToFront(m_gazecross);
			moveToFront(m_fixcross);
			m_isOpen = false;
			CEnvironment::Instance()->outputMessage("---Calibration Trial---");
		}
		else if (!m_isOpen && m_buttonPress) { // button confirmation received
			m_isOpen = true;
			m_gazecross->hide();
			m_fixcross->hide();
			endTrial();
			m_state = STATE_LOADING; // go to loading
		}
		else if (!m_isOpen) { // stabilize gazecross, wait for button press
			m_gazecross->pxSetPosition(X, Y);
		}

		break;
	case STATE_LOADING:
		if (m_isOpen) {
			cueCtr = 0; // reset landmark times
			cuePrf = 0;
			fixOn = 0;
			quit = 0;
			flashOn;
			rampOff = 0;
			saccOff = 0;
			saccOn = 0;
			stimOff = 0;

			m_isOpen = false;
			m_stimTimer.start(m_loadTime);

			// get parameter set and delete from list
			curr_trial = m_trial[m_trial.size() - 1];
			m_trial.pop_back();

			// calculate peripheral locations
			float fixXadeg = m_fixEcc * cos(curr_trial.fixAngle);
			float fixYadeg = m_fixEcc * sin(curr_trial.fixAngle);
			CConverter::Instance()->
				a2p(fixXadeg, fixYadeg, periphX, periphY);

			// reset saccade landing position
			landingX = NULL;
			landingY = NULL;

			// reset ramp frame counter
			currRampFrame = 0;
			
			//int instance = curr_trial.presTime + curr_trial.ecc + 2*curr_trial.presTime + 1;
			//int totalInstance = m_presTimeLevUsed * m_eccLevUsed;
			// if PEST is used, the contrast level is set here.
			if (!m_contrastLevUsed) {
				//curr_trial.contrast = m_pest[curr_trial.presTime][curr_trial.ecc]->getTestLvl();
				//curr_trial.contrast = m_pest[curr_trial.presTime]->getTestLvl();
				curr_trial.contrast = m_pest[curr_trial.instance]->getTestLvl();
				//curr_trial.contrast = m_contrast[0];
			}
			else { // otherwise pull randomly from the list
				curr_trial.contrast = m_contrast[rand() % m_contrastLevUsed];
			}

			// update the image + shader - contrast should range from 0 to .5
			
			m_noiseBackground->loadImage(curr_trial.background_img);
			m_imageShader->setUniformFloat("amplitude", 0); // reset the contrast of the grating
			m_imageShader->setUniformFloat("phase", curr_trial.phase); // reset the contrast of the grating
			m_imageShader->setUniformFloat("rcx", centerX); // central point of annulus
			m_imageShader->setUniformFloat("rcy", centerY);
			m_imageShader->setUniformFloat("eccentricity", m_ecc[curr_trial.ecc]);
			m_imageShader->setUniformFloat("SpatialFreq", m_spFreq[curr_trial.spFreq]);
			if (gratingType == "oriented grating") m_imageShader->setUniformFloat("tilt", curr_trial.present);

			CEnvironment::Instance()->outputMessage("==============================================");
			CEnvironment::Instance()->outputMessage("Trial %i - Remaining Trials - %i", m_currTrial, m_trial.size());
			CEnvironment::Instance()->outputMessage("background image: %s", curr_trial.background_img.c_str());
			CEnvironment::Instance()->outputMessage("contrast: %1.1f", curr_trial.contrast);
			CEnvironment::Instance()->outputMessage("Present %i, PresTime %1.1f, Ecc %1.1f", curr_trial.present, m_presTime[curr_trial.presTime],m_ecc[curr_trial.ecc]);
			CEnvironment::Instance()->outputMessage("SPATIAL FREQUENCY %1.1f CPD", m_spFreq[curr_trial.spFreq]);
		}
		else if (!m_isOpen && m_stimTimer.isExpired()) {
			m_isOpen = true;
			m_state = STATE_PERIPH_FIX;
		}
		break;

	case STATE_PERIPH_FIX:
		m_noiseBackground->show();
		if (m_isOpen) {
			startTrial();
			m_trialTimer.start(m_dropOutTime);
			cuePrf = m_trialTimer.getTime(); // cue start time

			if (!m_fixMode) { // if not fixation mode, show peripheral cue
				m_pCue->show();
				m_pCue->pxSetPosition(periphX, periphY);
				moveToFront(m_pCue);
			}
			m_crosshair->show();
			m_crosshair->pxSetPosition(centerX, centerY);
			moveToFront(m_crosshair);

			m_isOpen = false;
		}
		else if (m_fixMode || (getGazePxDist(periphX, periphY) < m_fixAreaPpx || debug)) {
			m_stimTimer.start(curr_trial.fixTime);
			fixOn = m_trialTimer.getTime();
			m_isOpen = true;
			m_state = STATE_SUSTAIN_FIX;
		}
		break;
	case STATE_SUSTAIN_FIX:
		m_noiseBackground->show();  // move this later
		if (!m_fixMode && // if gaze goes out of radius, go back to previous state
			(getGazePxDist(periphX, periphY) > m_fixAreaPpx && !debug) &&
			m_isOpen) {
			m_isOpen = false;
			m_state = STATE_PERIPH_FIX;
		}
		else if ((m_fixMode || (getGazePxDist(periphX, periphY) < m_fixAreaPpx || debug) ) &&
			m_stimTimer.isExpired() && // fixation time expired
			m_isOpen) {
			m_pCue->hide(); // hide peripheral cue;
			if (!m_fixMode) {
				m_cCue->hide(); // show central cue
				/*m_cCue->pxSetPosition(centerX, centerY);
				moveToFront(m_cCue);*/
			}

			cueCtr = m_trialTimer.getTime();
			m_isOpen = 0;
		}
		else if ((m_fixMode || (
			((CTriggers::all(Samples->triggers,
				Samples->samplesNumber, EOS_TRIG_1_EMEVENT) &&
			!CTriggers::any(Samples->triggers,
				Samples->samplesNumber, EOS_TRIG_1_BLINK)) ||
				(debug && rand() % 1000 == 0))
			)) &&
			!m_isOpen) {
			// SACCADE STARTED
			saccOn = m_trialTimer.getTime();
			m_stimTimer.start(m_saccDur);
			if (!m_fixMode) {
				m_cCue->hide();
			}

			if (m_flashed) {
				m_imageShader->setUniformFloat("amplitude", curr_trial.contrast*curr_trial.present);
				moveToFront(m_crosshair);
				flashOn = m_trialTimer.getTime();
			}
			m_isOpen = true;
			m_state = STATE_CENTRAL_FIX;
		}

		break;

	case STATE_CENTRAL_FIX:
		m_noiseBackground->show();  // move this later
		if (!m_stimTimer.isExpired()) {
			if (m_fixMode || // check for good central fixation
				(getGazePxDist(centerX, centerY) < m_fixAreaCpx &&
					!CTriggers::any(Samples->triggers,
						Samples->samplesNumber,
						EOS_TRIG_1_BADDATA) &&
					!CTriggers::any(
						Samples->triggers,
						Samples->samplesNumber,
						EOS_TRIG_1_EMEVENT) || debug)) {
				if (!m_flashed && currRampFrame < rampFrame) {
					currRampFrame = 0; // gaze has landed, start ramp
					m_isOpen = true;
				}
				else { // plateau or flash period, start full presentation
					m_stimTimer.start(m_presTime[curr_trial.presTime]);
					m_isOpen = false;
				}

				saccOff = m_trialTimer.getTime();
				moveToFront(m_crosshair);
				m_state = STATE_IMAGE;
			}
		}
		else { // trial timer ran out
			quit = m_trialTimer.getTime();
			m_trialTimer.stop();
			hideAllObjects();
			m_state = STATE_QUIT;
		}
		break;

	case STATE_IMAGE:
		m_noiseBackground->show();  // move this later
		if (m_stabilized) {
			if (landingX == NULL && landingY == NULL) {
				CStabilizer::Instance()->
					resetFilter(Samples->x1, Samples->y1);
				landingX = X;
				landingY = Y;
			}

			CStabilizer::Instance()->stabilize(Samples, Xstab, Ystab);
			xpos = Xstab - landingX + m_corrX;
			ypos = Ystab - landingY + m_corrY;

			m_crosshair->pxSetPosition(xpos, ypos);

			// change position of annulus in noise background
			m_imageShader->setUniformFloat("rcx", xpos);
			m_imageShader->setUniformFloat("rcy", ypos);
			storeTrialStream(0, xpos);
			storeTrialStream(1, ypos);
		}

		if (m_isOpen && currRampFrame < rampFrame) { // ramp phase
			currRampFrame++;
			float amp = curr_trial.contrast * currRampFrame / rampFrame * curr_trial.present;
			m_imageShader->setUniformFloat("amplitude", amp);
		}
		else if (m_isOpen && currRampFrame == rampFrame) { // transition from ramp to plateau
			m_imageShader->setUniformFloat("amplitude",curr_trial.contrast*curr_trial.present);
			rampOff = m_trialTimer.getTime();
			m_isOpen = false;
			m_stimTimer.start(m_presTime[curr_trial.presTime]);
		}
		else if (!m_isOpen && m_stimTimer.isExpired()) { // end plateau
			m_imageShader->setUniformFloat("amplitude", 0);
			m_crosshair->hide();
			stimOff = m_trialTimer.getTime();
			m_isOpen = true;
			m_state = STATE_WAIT_RESPONSE;
		}
		break;

	case STATE_WAIT_RESPONSE:
		m_noiseBackground->show();  // move this later
		if (m_isOpen) { // start response period
			//Beep(5000, 250);
			m_isOpen = false;
			m_rt.start(m_respTime);
		}
		else if (m_rt.isExpired() || m_buttonPress) {
			quit = m_trialTimer.getTime();
			m_isOpen = true;
			m_state = STATE_QUIT;
		}
		break;

	case STATE_QUIT:
		endTrial();

		
		hideAllObjects();
		if (curr_trial.rt >= 0) { // save data if trial did not time out
			saveData();
			m_trialTimer.stop();
		}
		else { // otherwise re-insert trial parameters into sequence and shuffle
			m_trial.push_back(curr_trial);
			random_shuffle(m_trial.begin(), m_trial.end());
		}

		// feedback
		if (m_feedback && curr_trial.resp != curr_trial.present) {
			Beep(2000, 150);
		}
		else {
			Beep(4000, 150);
		}

		// feed pest
		if (!m_contrastLevUsed && curr_trial.present && curr_trial.resp >= 0) {
			//m_pest[curr_trial.presTime][curr_trial.ecc]->addTrial(curr_trial.resp == curr_trial.present);
			//m_pest[curr_trial.presTime]->addTrial(curr_trial.resp == curr_trial.present);
			m_pest[curr_trial.instance]->addTrial(curr_trial.resp == curr_trial.present);
		}

		m_currTrial++;
		m_isOpen = true;



		if (!m_trial.empty()) {
			if (!m_fixMode && (m_currTrial % nRecal) == 0) {
				m_state = STATE_CALIBRATION;
			}
			else {
				m_state = STATE_LOADING;
			}
		}
		else {
			finalize();
		}

		m_buttonPress = false;
		break;

	}
}


///////////////////////////////////////////////////////////////////////////////////
void ExperimentBody::eventJoypad()
{

	switch (m_state)
	{
	case STATE_CALIBRATION:

		if (CDriver_Joypad::Instance()->getButtonStatus(CDriver_Joypad::JPAD_BUTTON_UP)) // moving the cursor up
		{
			m_corrY += Increment; //position of the cross
		}

		else if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_DOWN)) // moving the cursor down
		{
			m_corrY -= Increment;
		}

		else if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_RGHT)) // moving the cursor to the right
		{
			m_corrX += Increment;
		}

		else if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_LEFT)) // moving the cursor to the left
		{
			m_corrX -= Increment;
		}


		if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_R1))
		{ // finalize the response
			m_buttonPress = true; // click R1 for 
		}
		break;
	case STATE_WAIT_RESPONSE:
		if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_L1)) {
			curr_trial.resp = 0; // not present or left tilt
			curr_trial.rt = m_trialTimer.getTime();
		}

		if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_R1)) {
			curr_trial.resp = 1; // not present or left tilt
			curr_trial.rt = m_trialTimer.getTime();
		}

		if (curr_trial.resp >= 0) m_buttonPress = true;
		break;
	}

}

void ExperimentBody::eventKeyboard(unsigned char key, int x, int y) {
	switch (m_state)
	{
	case STATE_CALIBRATION:

		if (key == 'w' || key == 'W') // moving the cursor up
		{
			m_corrY += Increment; //position of the cross
		}

		else if (key == 's' || key == 'S') // moving the cursor down
		{
			m_corrY -= Increment;
		}

		else if (key == 'd' || key == 'D') // moving the cursor to the right
		{
			m_corrX += Increment;
		}

		else if (key == 'a' || key == 'A') // moving the cursor to the left
		{
			m_corrX -= Increment;
		}


		if (key == 'r' || key == 'R')
		{ // finalize the response
			m_buttonPress = true; // click R1 for 
		}
		break;
	case STATE_WAIT_RESPONSE:
		if (key == 'a' || key == 'A') {
			curr_trial.resp = 0; // not present or left tilt
			curr_trial.rt = m_trialTimer.getTime();
		}

		if (key == 'd' || key == 'D') {
			curr_trial.resp = 1; // not present or left tilt
			curr_trial.rt = m_trialTimer.getTime();
		}

		if (curr_trial.resp >= 0) m_buttonPress = true;
		break;
	}
}

void ExperimentBody::setParameters() {
	// determine the pixelAngle of the current settings
	pixelAngle = CConverter::Instance()->px2a(1.f); // arcmin per pixel
	CEnvironment::Instance()->outputMessage("initialize grating with pixelAngle = %1.3f", pixelAngle);

	//refreshRate = COGLEngine::Instance()->getFrameRate();

	// parameters for recalibration trials
	Increment = 10; // jump size of button press during recal trials (pixels)
	m_corrX = 0;
	m_corrY = 0;

	xpos = 0;
	ypos = 0;

	m_currTrial = 1;
	m_numTestCalibration = 0;
	// set TestCalibration = 1 so that the experiment will start with 
	// a recalibration trial
	m_numCompleted = 0;

	// General Settings	
	m_fixMode = m_paramsFile->getBoolean(CFG_FIX_MODE);
	m_stabilized = m_paramsFile->getBoolean(CFG_STABILIZED);
	m_feedback = m_paramsFile->getBoolean(CFG_FEEDBACK);
	m_pestStep = m_paramsFile->getFloat(CFG_PEST_STEP);
	m_pestTarget = m_paramsFile->getFloat(CFG_PEST_TARGET);
	m_flashed = m_paramsFile->getBoolean(CFG_FLASHED);

	m_sbj = m_paramsFile->getString(CFG_SUBJECT_NAME);
	m_backLum = m_paramsFile->getInteger(CFG_BACKGROUND_GRAY);
	m_contrastBack = m_paramsFile->getFloat(CFG_CONTRAST_BACK);

	// Timing
	m_rampTime = m_paramsFile->getInteger(CFG_RAMP_TIME);
	rampFrame = int(round(float(m_rampTime) / refreshRate));

	m_presTime[0] = m_paramsFile->getInteger(CFG_PRES_TIME_0);
	m_presTime[1] = m_paramsFile->getInteger(CFG_PRES_TIME_1);	

	m_respTime = m_paramsFile->getInteger(CFG_RESP_TIME);
	m_loadTime = m_paramsFile->getInteger(CFG_LOAD_TIME);
	m_fixTime = m_paramsFile->getInteger(CFG_FIX_TIME);
	m_fixTimeRng = m_paramsFile->getInteger(CFG_FIX_TIME_RNG);
	m_fixAngleOffset = m_paramsFile->getFloat(CFG_FIX_ANGLE_OFFSET);

	// Limiters
	m_presTimeLevUsed = m_paramsFile->getInteger(CFG_PRES_TIME_LEV);
	m_contrastLevUsed = m_paramsFile->getInteger(CFG_CONTRAST_LEV);
	m_fixLocNo = m_paramsFile->getInteger(CFG_FIX_LOC_NO);
	m_trialNo = m_paramsFile->getInteger(CFG_NUMB_TRIAL);
	m_imgNo = m_paramsFile->getInteger(CFG_NUMB_IMG);

	m_spatialFreqLev = m_paramsFile->getInteger(CFG_SPAT_FREQ_LEV);
	m_spFreq[0] = m_paramsFile->getInteger(CFG_SPAT_FREQ_0);
	m_spFreq[1] = m_paramsFile->getInteger(CFG_SPAT_FREQ_1);

	m_fixEcc = m_paramsFile->getFloat(CFG_FIX_ECC);

	// fixation areaP and C in arcmins
	m_fixAreaC = m_paramsFile->getFloat(CFG_FIX_AREA_C);
	m_fixAreaP = m_paramsFile->getFloat(CFG_FIX_AREA_P);

	m_contrast[0] = m_paramsFile->getFloat(CFG_CONTRAST_0);
	m_contrast[1] = m_paramsFile->getFloat(CFG_CONTRAST_1);
	m_contrast[2] = m_paramsFile->getFloat(CFG_CONTRAST_2);
	m_contrast[3] = m_paramsFile->getFloat(CFG_CONTRAST_3);
	m_contrast[4] = m_paramsFile->getFloat(CFG_CONTRAST_4);
	m_contrast[5] = m_paramsFile->getFloat(CFG_CONTRAST_5);
	m_contrast[6] = m_paramsFile->getFloat(CFG_CONTRAST_6);


	// Other Settings
	m_currTrial = 0;
	m_dropOutTime = 12000;
	m_corrX = m_corrY = 0;

	// chnage fixation areaP and C in pxls
	m_fixAreaCpx = CConverter::Instance()->
		a2px(m_fixAreaC);
	m_fixAreaPpx = CConverter::Instance()->
		a2px(m_fixAreaP);

	// Image Width and Height should be updated every time
	//	the size of the images used changes
	centerX = centerY = 0;


	// max time elapsing between saccade start and saccade end detection
	//	when this time is exceeded, trials gets aborted
	m_saccDur = 200; //800;//500;
	
	m_buttonPress = false;
	m_isOpen = true;

	awidth = m_paramsFile->getFloat(CFG_AWID);					// width of annulus
	gaussSD = m_paramsFile->getFloat(CFG_GSD);					// std of gaussian envelope
	gaussSD_2 = m_paramsFile->getFloat(CFG_GSD2);               // std of gaussian enevelope when ecc=0
	//debug = m_paramsFile->getInteger(CFG_DEBUG) == 1;				// debug mode?		
	debug = m_paramsFile->getBoolean(CFG_DEBUG);				// debug mode?	
	nRecal = m_paramsFile->getInteger(CFG_NRECAL);				// number of trials between recalibrations
	gratingType = m_paramsFile->getString(CFG_GT);				//raidal wave or oriented grating

	// to do: blocked eccentricity or truly pseudo-randomize?
	//eccentricity = m_paramsFile->getInteger(CFG_ECC);
	m_eccLevUsed = m_paramsFile->getInteger(CFG_ECC_LEV);
	m_ecc[0] = m_paramsFile->getInteger(CFG_ECC_0);
	m_ecc[1] = m_paramsFile->getInteger(CFG_ECC_1);
	m_ecc[2] = m_paramsFile->getInteger(CFG_ECC_2);	
	
	m_totalInstances = m_presTimeLevUsed * m_eccLevUsed*m_spatialFreqLev;
	//CEnvironment::Instance()->outputMessage("SPATIAL FREQUENCY %i CPD", spatialFreq);
	//CEnvironment::Instance()->outputMessage("ECCENTRICITY %i CPD", eccentricity); 
	CEnvironment::Instance()->outputMessage("Debug is %d", debug);
	CEnvironment::Instance()->outputMessage("FLASHED %i, RAMP TIME %i", m_flashed, m_rampTime);
}

void ExperimentBody::createObjects() {
	int fixationSize = 20;

	/* PEST initialization:
	if a pest.txt file is present in the data folder
	the pest is restored, otherwise a new one is started */
	if (!m_contrastLevUsed)
	{
		// TO DO: there should be a differerent file for all of the different conditions!!!!
		string pestFile = "./Data/" + m_sbj + "/pest.txt";
		ifstream ifile;
		ifile.open(pestFile, ios::in);
		if (ifile.good()) {
			ifile.close();
			//for (int i = 0; i < m_presTimeLevUsed; i++)
			for (int i = 0; i < m_totalInstances; i++) //TO DO:a nested loop for presLev and eccLevels
			{
				string pestFile = "./Data/" + m_sbj;
				//m_pest.push_back(vector<Pest*>());
				//m_pest[i].push_back(new Pest(pestFile.c_str()));
				m_pest.push_back(new Pest(pestFile.c_str()));
			}

			CEnvironment::Instance()->
				outputMessage("Loading PEST procedure started on: %s",
					Pest::getLabel());
		}
		else {
			// TO DO - make performance target a parameter?
			//for (int i = 0; i < m_presTimeLevUsed; i++) //TO DO:a nested loop for presLev and eccLevels
			//for (int i = 0; i < m_presTimeLevUsed; i++) //TO DO:a nested loop for presLev and eccLevels
			for (int i = 0; i < m_totalInstances; i++) //TO DO:a nested loop for presLev and eccLevels
			{
				//m_pest.push_back(vector<Pest*>());
				//m_pest[i].push_back(new Pest(m_contrast[0], 0.75, m_pestStep, 1, 1));
				m_pest.push_back(new Pest(m_contrast[0], m_pestTarget, m_pestStep, 1, 1));//*Nikunj* (startingLevel, targetP, stepSize, waldConstant, isLogarithmic)
			}

			//for (int i = 0; i < m_presTimeLevUsed; i++) //TO DO:a nested loop for presLev and eccLevels
			//{
			//	for (int j = 0; j < m_eccLevUsed; j++)
			//	{
			//		m_pest[i][j].push_back(vector<Pest*>());
			//		m_pest[i].push_back(new Pest(m_contrast[0], 0.75, m_pestStep, 1, 1));
			//		//m_pest.push_back(new Pest(m_contrast[i], 0.75, m_pestStep, 1, 0));//*Nikunj* (startingLevel, targetP, stepSize, waldConstant, isLogarithmic)
			//	}
			//}
		}
	}

	// create the fixation stimulus
	m_pCue = addObject(new CSolidPlane(190, 190, 190));
	m_pCue->pxSetSize(fixationSize, fixationSize);
	m_pCue->pxSetPosition(periphX, periphY);
	m_pCue->hide();

	m_cCue = addObject(new CSolidPlane(0, 0, 0));
	m_cCue->pxSetSize(fixationSize, fixationSize);
	m_cCue->pxSetPosition(centerX, centerY);
	m_cCue->hide();

	// boxes for the recalibration trials	
	// square versions of recalibration trial markers
	m_fixcross = addObject(new CSolidPlane(255, 255, 255));
	m_fixcross->pxSetSize(12, 12);
	m_gazecross = addObject(new CSolidPlane(0, 0, 0));
	m_gazecross->pxSetSize(12, 12);

	// cross hair - not sure if this will be feasible with the annulus
	m_crosshair = addObject(new CImagePlane("crosshair1.tga"));//*Nikunj* if changing it change it in the trial params as well so that the correct name of file can be stored.
	m_crosshair->enableTrasparency(true);
	m_crosshair->pxSetSize(200, 200);
	m_crosshair->hide();

	if (gratingType == "oriented grating") { // JI: 11-15-2018 - this option is not working yet
		m_imageShader = std::make_shared<CShader>("shaders/shader_main.vert", "shaders/shader_grating.frag");
		m_imageShader->setUniformFloat("tilt", 0);
	}
	else if (gratingType == "radial wave") {
		m_imageShader = std::make_shared<CShader>("shaders/shader_main.vert", "shaders/shader_radial_image.frag");
	}

	m_imageShader->setUniformFloat("gaussSD", gaussSD);
	m_imageShader->setUniformFloat("gaussSD2", gaussSD_2);
	m_imageShader->setUniformFloat("width", awidth);
	m_imageShader->setUniformFloat("pixelAngle", pixelAngle);
	//m_imageShader->setUniformFloat("eccentricity", eccentricity);
	//m_imageShader->setUniformFloat("SpatialFreq", spatialFreq);
	m_imageShader->setUniformFloat("amplitude", m_contrast[0]);
	m_imageShader->setUniformFloat("transparency", 1);
	m_imageShader->setUniformFloat("phase", 2); // will randomize later
	m_imageShader->setUniformFloat("imageScale", m_contrastBack); // will randomize later

	// Add shader to this object and enable it.
	// Enable trasparency for the shader.
	// apply a shader to a background noise image
	m_noiseBackground = addObject(new CImagePlane("noise/noise_1.bmp")); // load any file to start
	m_noiseBackground->pxSetPosition(0.f, 0.f);
	m_noiseBackground->show();

	m_noiseBackground->applyShader(m_imageShader);
	m_noiseBackground->enableShader(true);
	m_noiseBackground->enableTrasparency(true);
}



/******* createTrialSequnce: pre-randomize trials to pull from pseudorandomly *****/
void ExperimentBody::createTrialSequence() {
	// initialize single trials
	random_device rd;
	mt19937 mt(rd());
	uniform_int_distribution<int> catch_present(0, 1);
	uniform_int_distribution<int> eccentricity(0, m_eccLevUsed - 1);
	uniform_int_distribution<int> presentation(0, m_presTimeLevUsed - 1);
	uniform_int_distribution<int> spatialFreq(0, m_spatialFreqLev - 1);
	while (m_trial.size() < m_trialNo) {
		trialParams params;	

		params.seed = rd();
		// randomize stimulus parameters		
		params.present = catch_present(mt);;
		//params.present = rand() % 3 > 0; // doubles as tilt if oriented grating is used
		if (debug = 1)
			params.present = 1;		
		
		/*params.presLev = rand() % m_presTimeLevUsed;
		params.presTime = m_presTime[params.presLev];*/
		params.fixTime = m_fixTime - m_fixTimeRng / 2 + rand() % int(m_fixTimeRng);
		params.fixAngle = (m_fixAngleOffset - int(m_fixAngleOffset) +
			float((int(m_fixAngleOffset) + rand()) % m_fixLocNo)) *
			(PI * 2.0 / m_fixLocNo);
		params.phase = float(rand() % 100) / 100 * 3.14159;

		int background_idx = rand() % m_imgNo;
		params.background_img = "noise/noise_" + to_string(background_idx) + ".bmp";
		params.central_img = "crosshair1.tga";

		params.resp = -1; // response
		params.rt = -1; // reset reaction time

		params.presTime = presentation(mt);
		params.ecc = eccentricity(mt);
		params.spFreq = spatialFreq(mt);
		params.instance = (6 * params.presTime) + (2 * params.ecc) + (params.spFreq);

		//params.presTime = rand() % m_presTimeLevUsed; //*Nikunj* params.presTime stores the pres. level
		//params.ecc = rand() % m_eccLevUsed; //*Nikunj* params.ecc stores the ecc. level
		//params.spFreq = rand() % m_spatialFreqLev;//*Nikunj* params.spFreq stores the Spatial Frequency level
		/*params.instance = params.presTime + params.ecc + (2 * params.presTime);*/
		//params.instance = (6 * params.presTime) + (2 * params.ecc) + (params.spFreq);
		//params.totalInstances = m_presTimeLevUsed * m_eccLevUsed -1;	

		m_trial.push_back(params);
	}
}

/******* saveData: save trial variables and write them to some kind of output file ********/
void ExperimentBody::saveData() {

	CEnvironment::Instance()->outputMessage("Correct? %i", (curr_trial.resp == curr_trial.present));

	storeTrialVariable("Subject_Name", m_sbj);

	// monitor settings
	storeTrialVariable("backgroundLum", m_backLum);
	storeTrialVariable("screenw", m_pxWidth);
	storeTrialVariable("screenh", m_pxHeight);
	storeTrialVariable("screenr", refreshRate);
	storeTrialVariable("pixelAngle", pixelAngle);

	// experiment parameters
	storeTrialVariable("stabilized", m_stabilized);
	storeTrialVariable("flashed", m_flashed);
	storeTrialVariable("fixMode", m_fixMode);
	storeTrialVariable("rampTime", m_rampTime);
	storeTrialVariable("responseWaitTime", m_respTime);
	storeTrialVariable("fixTimeRange", m_fixTimeRng);
	storeTrialVariable("loadTime", m_loadTime);
	storeTrialVariable("dropOutTime", m_dropOutTime);
	storeTrialVariable("feedback", m_feedback);

	storeTrialVariable("pestStep", m_pestStep);
	storeTrialVariable("pestTarget", m_pestTarget);
	storeTrialVariable("saccDur", m_saccDur);
	storeTrialVariable("fixEcc", m_fixEcc);
	storeTrialVariable("fixAngleOffset", m_fixAngleOffset);
	storeTrialVariable("fixAreaC", m_fixAreaC);
	storeTrialVariable("fixAreaP", m_fixAreaP);
	storeTrialVariable("fixAreaCpx", m_fixAreaCpx);
	storeTrialVariable("fixAreaPpx", m_fixAreaPpx);

	// general stimulus parameters
	//storeTrialVariable("eccentricity", eccentricity); // to do: update
	storeTrialVariable("eccentricity", m_ecc[curr_trial.ecc]); //*Nikunj*
	storeTrialVariable("spatialFreq", m_spFreq[curr_trial.spFreq]);
	storeTrialVariable("randomDeviceSeed", curr_trial.seed);
	storeTrialVariable("pestInstance", curr_trial.instance);
	storeTrialVariable("gratingType", gratingType);
	storeTrialVariable("backgroundContrast", m_contrastBack);
	storeTrialVariable("annulusWidth", awidth);
	storeTrialVariable("annulusGaussSD", gaussSD);
	storeTrialVariable("gaborGaussSD", gaussSD_2);

	// screen locations
	storeTrialVariable("centerX", centerX);
	storeTrialVariable("centerY", centerY);
	storeTrialVariable("periphX", periphX);
	storeTrialVariable("periphY", periphY);

	// trial parameters
	storeTrialVariable("contrast", curr_trial.contrast);
	//storeTrialVariable("contrast", (curr_trial.contrast)/ (m_contrastBack*0.5));
	storeTrialVariable("present", curr_trial.present);
	storeTrialVariable("phase", curr_trial.phase);
	storeTrialVariable("fixAngle", curr_trial.fixAngle);
	storeTrialVariable("fixTime", curr_trial.fixTime);
	storeTrialVariable("presTime", m_presTime[curr_trial.presTime]);
	storeTrialVariable("trialid", m_currTrial);
	storeTrialVariable("backgroundImage", curr_trial.background_img);
	storeTrialVariable("centralImage", curr_trial.central_img);

	// trial responses
	storeTrialVariable("resp", curr_trial.resp);
	storeTrialVariable("responseTime", curr_trial.rt);
	storeTrialVariable("landingX", landingX);
	storeTrialVariable("landingY", landingY);

	storeTrialVariable("xoffset", m_corrX);
	storeTrialVariable("yoffset", m_corrY); //px
	
	// trial timing
	storeTrialVariable("cuePrf", cuePrf);
	storeTrialVariable("fixOn", fixOn);
	storeTrialVariable("cueCtr", cueCtr);
	storeTrialVariable("saccOn", saccOn);
	storeTrialVariable("saccOff", saccOff);
	storeTrialVariable("flashOn", flashOn);
	storeTrialVariable("rampOff", rampOff);
	storeTrialVariable("stimOff", stimOff);
	storeTrialVariable("quit", quit);

	m_numCompleted++;

	saveTrial("./Data/" + m_sbj);

	// keep track of the test calibration trials
	m_numTestCalibration++;
	// recalibration active at every nRecal trials
	if (m_numTestCalibration == nRecal)
	{
		m_numTestCalibration = 0;
	}

	if (m_trial.empty()) {
		declareFinished();
	}

	m_currTrial++;

	CEnvironment::Instance()->outputMessage("-----------------------------------------------------");
	m_state = STATE_LOADING;
}

///////////////////////////////////////////////////////////////////////////////////
void ExperimentBody::finalize()
{
	// TO DO: save PEST to file!!!
	string m_dFolder = "./Data/" + m_sbj;
	while (m_pest.size() != 0)
	{
		m_pest.front()->saveInstance(m_dFolder.c_str() );
		delete m_pest.front();
		m_pest.erase(m_pest.begin());
	}
}

float ExperimentBody::getGazePxDist(float pX, float pY) {
	float d = sqrt(pow(pX - X, 2) +
		pow(pY - Y, 2));
	return(d);
}
