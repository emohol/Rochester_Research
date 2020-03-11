#ifndef _MY_AUTO_CALIBRATOR_2_HPP
#define _MY_AUTO_CALIBRATOR_2_HPP

#include <emil/libc_common.hpp>
#include <emil/CTimer.hpp>
#include <emil/CDSPDataStream.hpp>
#include <emil/ExCalibrator.hpp>
#include <emil/CEOS.hpp>
#include <emil/CImagePlane.hpp>

/// Automatic calibration experiment
/** This experiment exists to calibrate the mapping between voltages as reported by the
eyetracker and arcminutes of visual angle. This is accomplished by displaying nine points
in a square grid onscreen, and instructing the subject to fixate on each one in turn
and capturing the gaze location in response to a button press when the subject feels
he/she is in the center of the fixation cross.  The experiment determines the fixation
location in voltage space, and computes regressions to determine the mapping into
angle space. */
class myAutoCalibration : public ExCalibrator
{
public:
	/// Default constructor
	/** @param pxWidth Width of the monitor in pixels
	@param pxHeight Height of the monitor in pixels
	@param hzRefreshRate Refresh rate of the monitor, in Hz */
	myAutoCalibration(int pxWidth, int pxHeight, int hzRefreshRate, bool monitor = false);

	/// Default destructor
	virtual ~myAutoCalibration() {};

	/// Standard event handlers

	/// @see CExperiment.finalize()
	virtual void finalize();

	/// @see CExperiment.eventKeyboard()
	virtual void eventKeyboard(unsigned char key, int x, int y);

	/// @see CExperiment.eventJoypad()
	virtual void eventJoypad();

	/// @see CExperiment.eventRender()
	virtual void eventRender(unsigned int frameCount, CEOSData* samples);

	/// @see CExperiment.initialize()
	virtual void initialize();

	void setPointOffset(unsigned int x, unsigned int y);

	void setBackgroundColor(unsigned char r, unsigned char g, unsigned char b)
	{
		m_backgroundColor.r = r;
		m_backgroundColor.g = g;
		m_backgroundColor.b = b;
	}

	void setBackgroundColor(EIS_RGB color)
	{
		m_backgroundColor = color;
	}


private:

	/// Calculate the fixation center
	/** To ensure that the true center is found, capture the gaze location upon subject's
	button press.*/
	bool calculateCenter(const vector<float>& dataVector, float threshold, float& center);

	/// Calculate array of quadrant indices 
	void calculateQuadrants(int centersIndices[]);

	/// Load an old calibration
	void loadOldCalibration();

	/// Initialize a new calibration
	void startNewCalibration();

	/// The current calibration state
	enum {
		STATE_QUESTION,
		STATE_BADLOAD,
		STATE_POINT,
		STATE_ERROR
	} m_state;

	/// The fixation dot
	CImagePlane* m_planeStimulus;

	/// File name of the stimulus used during the calibration
	std::string m_fileStimulus;

	/// The currently-displayed calibration point
	int m_currentPoint;

	/// Buffer to hold the raw eyetracker data
	std::vector<float> m_xData1;
	std::vector<float> m_yData1;
	std::vector<float> m_xData2;
	std::vector<float> m_yData2;
	/// A timer for pausing to allow subject to fixate
	CTimer m_timer;

	/// Pointer to the welcome banner
	CImagePlane* m_welcomeBanner;

	float m_paramCalThreshold; ///< Standard deviation threshold for the center-finding algorithm
	unsigned int m_paramMaxData; ///< Approximate maximum samples collected per point
	unsigned int m_paramPrepointWait; ///< Seconds to wait after displaying new point before collecting data
	unsigned int m_paramMinValidData;    ///< Minimum valid data to calculate center
	unsigned int m_paramReactionDelay;    ///<Offset into data for center calculation due to subject's reaction time
	unsigned int m_paramWindow; ///<Size of data-window to use in calibration

	unsigned int POINT_X_OFFSET_EXT = POINT_X_OFFSET;
	unsigned int POINT_Y_OFFSET_EXT = POINT_Y_OFFSET;

};

#endif    // _EX_AUTO_CALIBRATOR_2_HPP