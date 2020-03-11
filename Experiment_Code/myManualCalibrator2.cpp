#include "stdafx.h"
#include "myManualCalibrator2.h"
//#include "ljackuw.h"

#include <common/eos_common.hpp>

#include <emil/libc_common.hpp>
#include <emil/CEnvironment.hpp>
#include <emil/CEnvVariables.hpp>
#include <emil/CExceptions.hpp>
#include <emil/CEOS.hpp>
#include <emil/CFontEngine.hpp>
#include <emil/CMath.hpp>
#include <emil/CStabilizer.hpp>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Public methods

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
myManualCalibrator2::myManualCalibrator2(int pxWidth, int pxHeight, int hzRefreshRate, EOS_EYE eye, bool monitor) :
    ExCalibrator(pxWidth, pxHeight, hzRefreshRate, monitor)
{
    setExperimentName("Monocular Linear Manual Calibrator");

    // Retrieve the file name of the stimulus used to indicate
    // the eye position (red cross)
    m_fileStimulus = CEnvVariables::Instance()->getEYERISPathFile(CFG_N_CAL_EYEPOSITIONFILE);

    // File name of the cross stimulus used as reference marker
    m_fileCross = CEnvVariables::Instance()->getEYERISPathFile(CFG_N_CAL_STIMULUSFILE);
	m_fileCross = "images/whitecross.tga";
	m_fileStimulus = "images/redcross.tga";

	m_eye = eye;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myManualCalibrator2::eventRender(unsigned int, CEOSData* Samples)
{
    // If the controller is in digital mode, flash the screen white
    showAllObjects();

    if ((CDriver_Joypad::Instance()->getButtonStatus(CDriver_Joypad::JPAD_BUTTON_LEFT)) ||
		(m_key == 'a') || (m_key == 'A'))
    {
        m_xOffset -= m_xOffsetJog;
    }
	else if ((CDriver_Joypad::Instance()->getButtonStatus(CDriver_Joypad::JPAD_BUTTON_RGHT)) ||
		     (m_key == 'd') || (m_key == 'D'))
    {
        m_xOffset += m_xOffsetJog;
		m_key = 0;
    }

	if ((CDriver_Joypad::Instance()->getButtonStatus(CDriver_Joypad::JPAD_BUTTON_UP)) ||
		(m_key == 'w') || (m_key == 'W'))
    {
        m_yOffset += m_yOffsetJog;
    }
	else if ((CDriver_Joypad::Instance()->getButtonStatus(CDriver_Joypad::JPAD_BUTTON_DOWN)) ||
		(m_key == 's') || (m_key == 'S'))
    {
        m_yOffset -= m_yOffsetJog;
    }

    // Draw the offset and gain onscreen so we can watch them change
    glColor3f(255, 255, 255);
    print(CFontEngine::FONTS_ARIAL_14, static_cast<float>(-m_pxWidth / 2) + 20.0f, 
        static_cast<float>(-m_pxHeight / 2) + 70.0f,
		EOS_RENDER_MONITOR_FIRST, false,
        "X Offset: %.3f", m_xOffset);

    print(CFontEngine::FONTS_ARIAL_14, static_cast<float>(-m_pxWidth / 2) + 20.0f, 
        static_cast<float>(-m_pxHeight / 2) + 50.0f,
		EOS_RENDER_MONITOR_FIRST, false,
        "Y Offset: %.3f", m_yOffset);

    // Now, render the stimulus plane in the new, recalibrated location
    float x;
    float y;
    if (Samples->samplesNumber > 0)
    {
        // Stabilize the X
        CStabilizer::Instance()->stabilize(Samples, x, y, m_eye);

        // Check the limits
        x = (fabs(x) < (m_pxWidth / 2)) ? x : CMath::sign(x) * (m_pxWidth / 2);
        y = (fabs(y) < (m_pxHeight / 2)) ? y : CMath::sign(y) * (m_pxHeight / 2);

        print(
			CFontEngine::FONTS_ARIAL_BLACK_10, 
			0.0f, -(m_pxHeight / 2.0f) + 10, 
			EOS_RENDER_MONITOR_FIRST, 
			false, 
			"x %f y %f", x, y);

        m_planeStimulus->pxSetPosition(x, y);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myManualCalibrator2::eventKeyboard(unsigned char key, int, int)
{
	m_key = key;
	// If g is pressed then the keyboard status is reinitialized
	if (key == 'g' || key == 'G')
	{
		// Store the initial values
		m_xOffset = 0;
		m_yOffset = 0;
	}

	float *xA;
	float *yA;

	if (m_eye == EOS_EYE_1) {

		xA = m_xA1;
		yA = m_yA1;
	}
	else {
		xA = m_xA2;
		yA = m_yA2;
	}

	if ((key == 'w' || key == 'W') ||
		(key == 'a' || key == 'A') ||
		(key == 's' || key == 'S') ||
		(key == 'd' || key == 'D'))
	{
		*(xA + m_cursorPoint) = *(xA + m_cursorPoint) + m_xOffset;
		*(yA + m_cursorPoint) = *(yA + m_cursorPoint) + m_yOffset;

		m_xOffset = 0;
		m_yOffset = 0;

		switch (m_cursorPoint)
		{
		case 0:
			calculateMap(0, 1, m_eye);
			break;
		case 1:
			calculateMap(0, 2, m_eye);
			break;
		case 2:
			calculateMap(1, 2, m_eye);
			break;
		case 3:
			calculateMap(0, 1, m_eye);
			calculateMap(2, 3, m_eye);
			break;
		case 4:
			calculateMap(0, 4, m_eye);
			break;
		case 5:
			calculateMap(1, 2, m_eye);
			calculateMap(3, 4, m_eye);
			break;
		case 6:
			calculateMap(2, 3, m_eye);
			break;
		case 7:
			calculateMap(2, 4, m_eye);
			break;
		case 8:
			calculateMap(3, 4, m_eye);
			break;
		}

		calibrateDSP(m_eye);
	}

	if (key == 'r' || key == 'R')
	{
		m_cursorPoint++;
		
		if (m_cursorPoint >= POINTS_NUM)
		{
			m_cursorPoint = 0;
		}

		m_crossStimulus->pxSetPosition(m_xCoo[m_cursorPoint], m_yCoo[m_cursorPoint]);
		m_crossStimulus->enableTrasparency(true);


		m_xOffset = 0;
		m_yOffset = 0;

		saveCalibrationToFile();
	}
	else if (key == 'l' || key == 'L')
	{
		m_cursorPoint--;
		if (m_cursorPoint < 0)
		{
			m_cursorPoint = POINTS_NUM;
		}

		m_crossStimulus->pxSetPosition(m_xCoo[m_cursorPoint], m_yCoo[m_cursorPoint]);
		m_crossStimulus->enableTrasparency(true);

		m_xOffset = 0;
		m_yOffset = 0;

		saveCalibrationToFile();
	}

	if (key == 'f' || key == 'F')
	{
		saveCalibrationToFile();
		declareFinished();
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myManualCalibrator2::eventJoypad()
{
    // If SELECT is pressed then the joypad status is reinitialized
    if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_SLCT))
    {
        // Store the initial values
        m_xOffset = 0;
        m_yOffset = 0;
    }

	float *xA;
	float *yA;

	if (m_eye == EOS_EYE_1) {

		xA = m_xA1;
		yA = m_yA1;
	} else {
		xA = m_xA2;
		yA = m_yA2;
	}

    if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_LEFT) ||
        CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_RGHT) ||
        CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_UP) ||
        CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_DOWN))
    {
		*(xA + m_cursorPoint) = *(xA + m_cursorPoint) + m_xOffset;
		*(yA + m_cursorPoint) = *(yA + m_cursorPoint) + m_yOffset;

        m_xOffset = 0;
        m_yOffset = 0;

        switch (m_cursorPoint)
        {
        case 0:
            calculateMap(0, 1, m_eye);
            break;
        case 1:
			calculateMap(0, 2, m_eye);
            break;
        case 2:
			calculateMap(1, 2, m_eye);
            break;
        case 3:
			calculateMap(0, 1, m_eye);
			calculateMap(2, 3, m_eye);
            break;
        case 4:
			calculateMap(0, 4, m_eye);
            break;
        case 5:
			calculateMap(1, 2, m_eye);
			calculateMap(3, 4, m_eye);
            break;
        case 6:
			calculateMap(2, 3, m_eye);
            break;
        case 7:
			calculateMap(2, 4, m_eye);
            break;
        case 8:
			calculateMap(3, 4, m_eye);
            break;
        }

        calibrateDSP(m_eye);
    } 
	
	if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_R1))
    {
        m_cursorPoint++;

        if (m_cursorPoint >= POINTS_NUM)
        {
            m_cursorPoint = 0;
        }

        m_crossStimulus->pxSetPosition(m_xCoo[m_cursorPoint], m_yCoo[m_cursorPoint]);
        m_crossStimulus->enableTrasparency(true);

        m_xOffset = 0;
        m_yOffset = 0;

        saveCalibrationToFile();
    }
    else if (CDriver_Joypad::Instance()->getButtonReleased(CDriver_Joypad::JPAD_BUTTON_L1))
    {
        m_cursorPoint--;
        if (m_cursorPoint < 0)
        {
            m_cursorPoint = POINTS_NUM;
        }

        m_crossStimulus->pxSetPosition(m_xCoo[m_cursorPoint], m_yCoo[m_cursorPoint]);
        m_crossStimulus->enableTrasparency(true);

        m_xOffset = 0;
        m_yOffset = 0;

        saveCalibrationToFile();
    }

    if (CDriver_Joypad::Instance()->getButtonPressed(CDriver_Joypad::JPAD_BUTTON_STRT))
    {
        saveCalibrationToFile();
        declareFinished();
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myManualCalibrator2::finalize()
{
    endTrial();

    calibrateDSP();
    saveCalibrationToFile();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void myManualCalibrator2::initialize()
{
    // Call its ancestor
    ExCalibrator::initialize();

    m_cursorPoint = 0;

    // Set the initial point coordinates
    for (int i = 0; i < POINTS_NUM; i++)
    {
        m_xCoo[i] = static_cast<short>(((i % 3) - 1) * POINT_X_OFFSET_EXT);
        m_yCoo[i] = static_cast<short>(-((i / 3) - 1) * POINT_Y_OFFSET_EXT);
    }

    for (unsigned int i = 0; i < POINTS_NUM; i++)
    {
        m_xA1[i] = CConverter::Instance()->px2a(m_xCoo[i]);
        m_yA1[i] = CConverter::Instance()->py2a(m_yCoo[i]);

		m_xA2[i] = CConverter::Instance()->px2a(m_xCoo[i]);
		m_yA2[i] = CConverter::Instance()->py2a(m_yCoo[i]);
    }

    // Set the eye position stimulus (red cross)
	if ((m_eye == EOS_EYE_2) & isBinocular() & isDualMonitor()) {
		m_planeStimulus = addObject(new CImagePlane(m_fileStimulus));
		m_planeStimulus->pxSetSize(20, 20);
		m_planeStimulus->setMonitor(EOS_RENDER_MONITOR_SECOND);
	}
	else 
	{
		m_planeStimulus = addObject(new CImagePlane(m_fileStimulus));
		m_planeStimulus->pxSetSize(20, 20);
		m_planeStimulus->setMonitor(EOS_RENDER_MONITOR_FIRST);
	}
    m_planeStimulus->pxSetPosition(0, 0);
    m_planeStimulus->enableTrasparency(true);

    // Set the initial point coordinates
	if ((m_eye == EOS_EYE_2) & isBinocular() & isDualMonitor()) {
		m_crossStimulus = addObject(new CImagePlane(m_fileCross));
		m_crossStimulus->pxSetSize(20, 20);
		m_crossStimulus->setMonitor(EOS_RENDER_MONITOR_SECOND);
	}
	else 
	{
		m_crossStimulus = addObject(new CImagePlane(m_fileCross));
		m_crossStimulus->pxSetSize(20, 20);
		m_crossStimulus->setMonitor(EOS_RENDER_MONITOR_FIRST);
	}
    m_crossStimulus->pxSetPosition(m_xCoo[m_cursorPoint], m_yCoo[m_cursorPoint]);
    m_crossStimulus->enableTrasparency(true);


	moveToFront(m_planeStimulus);
    loadCalibrationFromFile();

    // Initial values
    m_xOffset = 0;
    m_yOffset = 0;

    // Set the jog values
    float OffsetJog = CEnvVariables::Instance()->getFloat(CFG_N_CAL_OFFSETJOG);

    m_xOffsetJog = OffsetJog/m_refreshRate;
    m_yOffsetJog = OffsetJog/m_refreshRate;

    disable(CExperiment::EIS_JP_STRT);
    disable(CExperiment::EIS_PHOTOCELL);
    enable(CExperiment::EIS_STAT1);

    // Prepare the background
    COGLEngine::Instance()->setBackgroundColor(127, 127, 127);

    // Start the calibrated trial without recording
    startTrial(false, true);
}

void myManualCalibrator2::setPointOffset(unsigned int x, unsigned int y)
{
	POINT_X_OFFSET_EXT = x;
	POINT_Y_OFFSET_EXT = y;
}
