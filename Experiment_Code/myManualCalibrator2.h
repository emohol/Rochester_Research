#ifndef _MY_MANUAL_CALIBRATOR2_HPP
#define _MY_MANUAL_CALIBRATOR2_HPP

#include <emil/libc_common.hpp>
#include <emil/CExperiment.hpp>
#include <emil/CDSPDataStream.hpp>
#include <emil/ExCalibrator.hpp>
#include <emil/CImagePlane.hpp>

/// Manual calibration experiment
/** This experiment exists to further calibrate the mapping between voltages as reported by the
    eyetracker and arcminutes of visual angle. This is accomplished by letting the subject
    manually change gain and offset of the movements of a stabilized red X on few default
    points. */
class myManualCalibrator2 : public ExCalibrator
{
public:
    /// Default constructor
    /** @param pxWidth Width of the monitor in pixels
        @param pxHeight Height of the monitor in pixels
        @param hzRefreshRate Refresh rate of the monitor, in Hz */
    myManualCalibrator2(int pxWidth, int pxHeight, int hzRefreshRate, EOS_EYE eye, bool dualmonitor=false);

    /// Default destructor
    virtual ~myManualCalibrator2() {};

    /// Standard event handlers

    /// @see CExperiment.finalize()
    virtual void finalize();

    /// @see CExperiment.eventJoypad()
    virtual void eventJoypad();

	virtual void eventKeyboard(unsigned char key, int, int);

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


    // Some constants
    enum {
        LH_EYE_MIN = 90,
        LH_EYE_MAX = 160,

        LV_EYE_MIN = 90,
        LV_EYE_MAX = 160,

        RH_EYE_MIN = 90,
        RH_EYE_MAX = 160,

        RV_EYE_MIN = 90,
        RV_EYE_MAX = 160
    };

protected:
    /// The fixation dot
    CImagePlane* m_planeStimulus;

    ///The fixation cross
    CImagePlane* m_crossStimulus;

    /// File name of the stimulus used during the calibration
    std::string m_fileStimulus;

    /// File name of the cross stimulus used as reference marker
    std::string m_fileCross;

    int m_cursorPoint;

    float m_xOffset;
    float m_yOffset;

    float m_xOffsetJog;
    float m_yOffsetJog;

	EOS_EYE m_eye = EOS_EYE_1;

	char m_key=0;

	unsigned int POINT_X_OFFSET_EXT = POINT_X_OFFSET;
	unsigned int POINT_Y_OFFSET_EXT = POINT_Y_OFFSET;
};

#endif  // _EX_MANUAL_CALIBRATOR2_HPP
