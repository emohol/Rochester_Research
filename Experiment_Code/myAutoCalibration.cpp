#include "stdafx.h"
#include "myAutoCalibration.h"
#include "ljackuw.h"

#include <common/eos_common.hpp>

#include <emil/libc_common.hpp>
#include <emil/CMath.hpp>
#include <emil/CEnvironment.hpp>
#include <emil/CExceptions.hpp>
#include <emil/CEOS.hpp>
#include <emil/CEnvVariables.hpp>
#include <emil/CFontEngine.hpp>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Private methods

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool myAutoCalibration::calculateCenter(const vector<float>& dataVector, float threshold, float& center)
{
	// Return immediately if the input vector is just too small
	if (static_cast<int>(dataVector.size()) < m_paramMinValidData + m_paramReactionDelay)
	{
		CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_ERROR, "Vector too small");
		return false;
	}

	vector<float>::const_iterator begin = dataVector.end();
	vector<float>::const_iterator end = dataVector.end();

	begin -= m_paramReactionDelay + m_paramMinValidData;
	end -= m_paramReactionDelay;

	double stdDev = CMath::stddev(begin, end);
	if (stdDev < threshold)
	{
		center = CMath::mean(begin, end);
		return true;
	}

	CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_ERROR,
		"Unable to meet standard deviation threshold");

	return false;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::loadOldCalibration()
{
	if (loadCalibrationFromFile())
	{
		calibrateDSP();
		saveCalibrationToFile();
		declareFinished();
	}
	else
	{
		CDriver_Joypad::Instance()->clearButtons();
		m_state = STATE_BADLOAD;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::startNewCalibration()
{
	m_welcomeBanner->hide();
	COGLEngine::Instance()->setBackgroundColor(m_backgroundColor);

	CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_HIGHLIGHT, "Automatic calibration started");

	m_currentPoint = 0;
	m_state = STATE_POINT;

	// Activate the sample display
	enable(CExperiment::EIS_STAT1);

	// Show the first cross
	m_planeStimulus->pxSetPosition(m_xCoo[m_currentPoint], m_yCoo[m_currentPoint]);
	m_planeStimulus->enableTrasparency(true);
	m_planeStimulus->show();

	// Clear all previous calibration data
	m_xData1.clear();
	m_xData1.reserve(m_paramMaxData + 25);
	m_yData1.clear();
	m_yData1.reserve(m_paramMaxData + 25);


	m_xCenters1.clear();
	m_xCenters1.resize(POINTS_NUM);
	m_yCenters1.clear();
	m_yCenters1.resize(POINTS_NUM);

	bool binocular = CEOS::Instance()->isEOSBinocular();

	if (binocular) {
		m_xCenters2.clear();
		m_xCenters2.resize(POINTS_NUM);
		m_yCenters2.clear();
		m_yCenters2.resize(POINTS_NUM);

		m_xData2.clear();
		m_xData2.reserve(m_paramMaxData + 25);
		m_yData2.clear();
		m_yData2.reserve(m_paramMaxData + 25);
	}
	// Calculate the visual angles of points
	for (unsigned int i = 0; i < POINTS_NUM; i++)
	{
		m_xA1[i] = CConverter::Instance()->px2a(m_xCoo[i]);
		m_yA1[i] = CConverter::Instance()->py2a(m_yCoo[i]);
	}

	if (binocular) {
		for (unsigned int i = 0; i < POINTS_NUM; i++)
		{
			m_xA2[i] = CConverter::Instance()->px2a(m_xCoo[i]);
			m_yA2[i] = CConverter::Instance()->py2a(m_yCoo[i]);
		}
	}

	// Start data recording
	startTrial();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public methods

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
myAutoCalibration::myAutoCalibration(int pxWidth, int pxHeight, int hzRefreshRate, bool monitor) :
ExCalibrator(pxWidth, pxHeight, hzRefreshRate, monitor)
{
	setExperimentName("Monocular Linear Auto Calibrator");

	// Retrieve the name of the file that contains the stimulus to visualize
	// during the calibration
	m_fileStimulus = CEnvVariables::Instance()->getEYERISPathFile(CFG_N_CAL_STIMULUSFILE);
	m_fileStimulus = "images/whitecross.tga";

	// Initialize all parameters
	m_paramCalThreshold = CEnvVariables::Instance()->getFloat(CFG_N_CAL_CALSTDDEV);
	m_paramMaxData = CEnvVariables::Instance()->getInteger(CFG_N_CAL_MAXDATA);
	m_paramPrepointWait = CEnvVariables::Instance()->getInteger(CFG_N_CAL_PREPOINTWAIT);
	m_paramMinValidData = CEnvVariables::Instance()->getInteger(CFG_N_CAL_MINVALIDDATA);
	m_paramReactionDelay = CEnvVariables::Instance()->getInteger(CFG_N_CAL_REACTIONDELAY);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::eventKeyboard(unsigned char key, int, int)
{
	if (m_state == STATE_QUESTION)
	{
		// Wait for an answer
		if (key == 'l' || key == 'L')
		{
			loadOldCalibration();
		}
		else if (key == 'r' || key == 'R')
		{
			startNewCalibration();
		}
	}
	else if (m_state == STATE_POINT)
	{
		if (key == 'l' || key == 'L')
		{
			if (m_currentPoint != 0)
			{
				m_currentPoint--;

				// Prepare to collect the data
				m_xData1.clear();
				m_xData1.reserve(m_paramMaxData + 25);
				m_yData1.clear();
				m_yData1.reserve(m_paramMaxData + 25);

				if (isBinocular()) {
					m_xData2.clear();
					m_xData2.reserve(m_paramMaxData + 25);
					m_yData2.clear();
					m_yData2.reserve(m_paramMaxData + 25);
				}

				// Show the previous cross
				m_planeStimulus->pxSetPosition(m_xCoo[m_currentPoint], m_yCoo[m_currentPoint]);
				m_planeStimulus->show();
			}
		}
		else if (key == 'r' || key == 'R')
		{
			// End recording
			endTrial();
			bool error = false;
			bool binocular = CEOS::Instance()->isEOSBinocular();

			// Try to calculate center
			if (!calculateCenter(m_xData1, m_paramCalThreshold, m_xCenters1[m_currentPoint]) ||
				!calculateCenter(m_yData1, m_paramCalThreshold, m_yCenters1[m_currentPoint]))
			{
				CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM,
					"Calibration in point %d failed at EYE 1", m_currentPoint + 1);

				COGLEngine::Instance()->setBackgroundColor(255, 0, 0);
				m_timer.start(1000.f);
				m_state = STATE_ERROR;
				error = true;
			}

			if (binocular) {
				if (!calculateCenter(m_xData2, m_paramCalThreshold, m_xCenters2[m_currentPoint]) ||
					!calculateCenter(m_yData2, m_paramCalThreshold, m_yCenters2[m_currentPoint]))
				{
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM,
						"Calibration in point %d failed at EYE 2", m_currentPoint + 1);

					COGLEngine::Instance()->setBackgroundColor(255, 0, 0);
					m_timer.start(1000.f);
					m_state = STATE_ERROR;
					error = true;
				}
			}

			if (!error)
			{
				// Print the results
				if (binocular) {
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "     Point %d: center (%f, %f) (%f, %f)",
						m_currentPoint + 1,
						m_xCenters1[m_currentPoint], m_yCenters1[m_currentPoint],
						m_xCenters2[m_currentPoint], m_yCenters2[m_currentPoint]);
				}
				else {
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "     Point %d: center (%f, %f)",
						m_currentPoint + 1, m_xCenters1[m_currentPoint], m_yCenters1[m_currentPoint]);
				}
				// Write the data out to file if enabled
				if (CEnvVariables::Instance()->getBoolean(CFG_N_CAL_SAVEPOINTS))
				{
					saveTrial(m_dirOutput, "AUTOCAL");
				}


				// If all points have been calibrated, calibrate the DSP and exit
				if (m_currentPoint == POINTS_NUM)
				{
					// Calculate the interpolation map
					calculateMap(0, 4);
					if (binocular) calculateMap(0, 4, EOS_EYE_2);

					calibrateDSP();
					saveCalibrationToFile();
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "");
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_HIGHLIGHT, "Calibration successful");
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "");

					declareFinished();
				}
				else
				{
					// Otherwise, go to the next calibration point
					// Prepare to collect the data
					m_xData1.clear();
					m_xData1.reserve(m_paramMaxData + 25);
					m_yData1.clear();
					m_yData1.reserve(m_paramMaxData + 25);

					if (isBinocular()) {
						m_xData2.clear();
						m_xData2.reserve(m_paramMaxData + 25);
						m_yData2.clear();
						m_yData2.reserve(m_paramMaxData + 25);
					}

					// Show the next cross
					m_currentPoint++;
					m_planeStimulus->pxSetPosition(m_xCoo[m_currentPoint], m_yCoo[m_currentPoint]);
					m_planeStimulus->show();

					// Start data recording
					startTrial();
				}
			}
		}
	}
	else if (m_state == STATE_BADLOAD)
	{
		m_state = STATE_QUESTION;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::eventJoypad()
{
	if (m_state == STATE_QUESTION)
	{
		if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_L1))
		{
			loadOldCalibration();
		}
		else
		{
			if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_R1))
			{
				startNewCalibration();
			}
		}
	}
	else if (m_state == STATE_POINT)
	{
		if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_L1))
		{
			CEnvironment::Instance()->outputMessage(
				"L1 button released");
			if (m_currentPoint != 0)
			{
				m_currentPoint--;

				// Prepare to collect the data
				m_xData1.clear();
				m_xData1.reserve(m_paramMaxData + 25);
				m_yData1.clear();
				m_yData1.reserve(m_paramMaxData + 25);

				if (isBinocular()) {
					m_xData2.clear();
					m_xData2.reserve(m_paramMaxData + 25);
					m_yData2.clear();
					m_yData2.reserve(m_paramMaxData + 25);
				}

				// Show the previous cross
				m_planeStimulus->pxSetPosition(m_xCoo[m_currentPoint], m_yCoo[m_currentPoint]);
				m_planeStimulus->show();
			}
		}
		else if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_R1))
		{
			// End recording
			endTrial();
			CEnvironment::Instance()->outputMessage(
				"R1 button released");
			bool error = false;
			bool binocular = CEOS::Instance()->isEOSBinocular();

			// Try to calculate center
			if (!calculateCenter(m_xData1, m_paramCalThreshold, m_xCenters1[m_currentPoint]) ||
				!calculateCenter(m_yData1, m_paramCalThreshold, m_yCenters1[m_currentPoint]))
			{
				CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM,
					"Calibration in point %d failed at EYE 1", m_currentPoint + 1);

				COGLEngine::Instance()->setBackgroundColor(255, 0, 0);
				m_timer.start(1000.f);
				m_state = STATE_ERROR;
				error = true;
			}

			if (binocular) {
				if (!calculateCenter(m_xData2, m_paramCalThreshold, m_xCenters2[m_currentPoint]) ||
					!calculateCenter(m_yData2, m_paramCalThreshold, m_yCenters2[m_currentPoint]))
				{
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM,
						"Calibration in point %d failed at EYE 2", m_currentPoint + 1);

					COGLEngine::Instance()->setBackgroundColor(255, 0, 0);
					m_timer.start(1000.f);
					m_state = STATE_ERROR;
					error = true;
				}
			}

			if (!error)
			{
				// Print the results
				if (binocular) {
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "     Point %d: center (%f, %f) (%f, %f)",
						m_currentPoint + 1,
						m_xCenters1[m_currentPoint], m_yCenters1[m_currentPoint],
						m_xCenters2[m_currentPoint], m_yCenters2[m_currentPoint]);
				}
				else {
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "     Point %d: center (%f, %f)",
						m_currentPoint + 1, m_xCenters1[m_currentPoint], m_yCenters1[m_currentPoint]);
				}
				// Write the data out to file if enabled
				if (CEnvVariables::Instance()->getBoolean(CFG_N_CAL_SAVEPOINTS))
				{
					saveTrial(m_dirOutput, "AUTOCAL");
				}


				// If all points have been calibrated, calibrate the DSP and exit
				if (m_currentPoint == POINTS_NUM-1)
				{
					// Calculate the interpolation map
					calculateMap(0, 4);
					if (binocular) calculateMap(0, 4, EOS_EYE_2);

					calibrateDSP();
					saveCalibrationToFile();
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "");
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_HIGHLIGHT, "Calibration successful");
					CEnvironment::Instance()->outputMessage(CEnvironment::ENV_MSG_SYSTEM, "");

					declareFinished();
				}
				else
				{
					// Otherwise, go to the next calibration point
					// Prepare to collect the data
					m_xData1.clear();
					m_xData1.reserve(m_paramMaxData + 25);
					m_yData1.clear();
					m_yData1.reserve(m_paramMaxData + 25);

					if (isBinocular()) {
						m_xData2.clear();
						m_xData2.reserve(m_paramMaxData + 25);
						m_yData2.clear();
						m_yData2.reserve(m_paramMaxData + 25);
					}

					// Show the next cross
					m_currentPoint++;
					m_planeStimulus->pxSetPosition(m_xCoo[m_currentPoint], m_yCoo[m_currentPoint]);
					m_planeStimulus->show();

					// Start data recording
					startTrial();
					CEnvironment::Instance()->outputMessage(
						"trial started");
				}
			}
		}
	}
	else
	{
		if (m_state == STATE_BADLOAD && CDriver_Joypad::Instance()->anyButtonsReleased())
		{
			m_state = STATE_QUESTION;
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::eventRender(unsigned int, CEOSData* samples)
{
	switch (m_state)
	{
	case STATE_QUESTION: // Show the first panel

		m_welcomeBanner->renderObject();

		glColor3d(255, 255, 255);
		printCentered(CFontEngine::FONTS_ARIAL_BLACK_16, 0.0f, 50.0f,
			EOS_RENDER_MONITOR_BOTH, false,
			"Automatic calibration procedure");
		printCentered(CFontEngine::FONTS_ARIAL_14, 0.0f, -55.0f,
			EOS_RENDER_MONITOR_BOTH, false,
			"Press \"L\" or left triggers to load the old calibration from file");
		printCentered(CFontEngine::FONTS_ARIAL_14, 0.0f, -80.0f,
			EOS_RENDER_MONITOR_BOTH, false,
			"Press \"R\" or right triggers to start a new calibration");
		break;

	case STATE_BADLOAD: // The calibration data load was unsuccessful

		m_welcomeBanner->renderObject();

		glColor3d(255, 0, 0);
		printCentered(CFontEngine::FONTS_ARIAL_14, 0.0f, -55.0f, EOS_RENDER_MONITOR_BOTH, false, "Unable to load previous calibration");

		glColor3d(255, 255, 255);
		printCentered(CFontEngine::FONTS_ARIAL_14, 0.0f, -80.0f, EOS_RENDER_MONITOR_BOTH, false, "Press any key or button to continue");
		break;

	case STATE_ERROR:

		if (m_timer.isExpired())
		{
			COGLEngine::Instance()->setBackgroundColor(m_backgroundColor);

			m_xData1.clear();
			m_xData1.reserve(m_paramMaxData + 25);
			m_yData1.clear();
			m_yData1.reserve(m_paramMaxData + 25);
			
			if (isBinocular()) {
				m_xData2.clear();
				m_xData2.reserve(m_paramMaxData + 25);
				m_yData2.clear();
				m_yData2.reserve(m_paramMaxData + 25);
			}

			m_state = STATE_POINT;
			startTrial();
		}
		break;

	case STATE_POINT: // If we are recording data at the point

		// Add collected samples to our buffer
		for (int i = 0; i < samples->samplesNumber; i++)
		{
			m_xData1.push_back(samples->vchannels[i][EOS_VCHANNEL_X1]);
			m_yData1.push_back(samples->vchannels[i][EOS_VCHANNEL_Y1]);
			m_xData2.push_back(samples->vchannels[i][EOS_VCHANNEL_X2]);
			m_yData2.push_back(samples->vchannels[i][EOS_VCHANNEL_Y2]);
		}
		break;
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::finalize()
{
	delete m_welcomeBanner;
	// Call its ancestor
	ExCalibrator::finalize();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myAutoCalibration::initialize()
{
	// Call its ancestor
	ExCalibrator::initialize();

	// Initialize the color background
	COGLEngine::Instance()->setBackgroundColor(0, 0, 50);

	// Create the stimulus plane, hidden for now
	m_planeStimulus = addObject(new CImagePlane(m_fileStimulus));
	m_planeStimulus->pxSetSize(20, 20);
	m_planeStimulus->setMonitor(EOS_RENDER_MONITOR_BOTH);
	m_planeStimulus->hide();

	// Set the initial point coordinates
	for (int i = 0; i < POINTS_NUM; i++)
	{
		m_xCoo[i] = static_cast<short>(((i % 3) - 1) * POINT_X_OFFSET_EXT);
		m_yCoo[i] = static_cast<short>(-((i / 3) - 1) * POINT_Y_OFFSET_EXT);
	}

	// Set the initial state
	m_state = STATE_QUESTION;

	// Load the welcome banner
	m_welcomeBanner = new CImagePlane(CEnvVariables::Instance()->getDirectory(CFG_N_GENERAL_ROOTDIR) + "system/images/win-frame.tga");
	m_welcomeBanner->setMonitor(EOS_RENDER_MONITOR_BOTH);
	m_welcomeBanner->pxSetPosition(0, 0);
	m_welcomeBanner->show();

}

void myAutoCalibration::setPointOffset(unsigned int x, unsigned int y)
{
	POINT_X_OFFSET_EXT = x;
	POINT_Y_OFFSET_EXT = y;
}
