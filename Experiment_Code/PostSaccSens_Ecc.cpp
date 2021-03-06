
// PostSaccSens_Ecc.cpp : Defines the class behaviors for the application.
//

#include "stdafx.h"
#include "PostSaccSens_Ecc.h"
#include "ExperimentBody.h"

#include "myAutoCalibration.h"
#include "myManualCalibrator2.h"

#include <emil-console\CEMILConsole.hpp>
#include <emil-console\CWinEMIL.hpp>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CEyeChartEApp

BEGIN_MESSAGE_MAP(CPostSaccSens_EccApp, CWinApp)
	ON_COMMAND(ID_HELP, &CWinApp::OnHelp)
END_MESSAGE_MAP()


// CEyeChartEApp construction

CPostSaccSens_EccApp::CPostSaccSens_EccApp()
{
	// support Restart Manager
	m_dwRestartManagerSupportFlags = AFX_RESTART_MANAGER_SUPPORT_RESTART;

	// General Settings
	m_paramsFile.addVariable(CFG_FIX_MODE, true);
	m_paramsFile.addVariable(CFG_STABILIZED, false);
	m_paramsFile.addVariable(CFG_FLASHED, false);
	m_paramsFile.addVariable(CFG_FEEDBACK, false);
	m_paramsFile.addVariable(CFG_PEST_STEP, 1);
	m_paramsFile.addVariable(CFG_PEST_TARGET, 1);

	m_paramsFile.addVariable(CFG_SUBJECT_NAME, std::string(""));
	
	m_paramsFile.addVariable(CFG_SUBJECT_NAME, std::string(""));
	m_paramsFile.addVariable(CFG_BACKGROUND_GRAY, 127);
	m_paramsFile.addVariable(CFG_SCREEN_W, "1920");
	m_paramsFile.addVariable(CFG_SCREEN_H, "1440");
	m_paramsFile.addVariable(CFG_SCREEN_R, "60");

	// Timing 
	m_paramsFile.addVariable(CFG_RAMP_TIME, 5000);

	m_paramsFile.addVariable(CFG_PRES_TIME_0, 5000);
	m_paramsFile.addVariable(CFG_PRES_TIME_1, 10000);
	m_paramsFile.addVariable(CFG_PRES_TIME_2, 10000);
	m_paramsFile.addVariable(CFG_PRES_TIME_3, 10000);
	m_paramsFile.addVariable(CFG_PRES_TIME_4, 10000);
	m_paramsFile.addVariable(CFG_PRES_TIME_5, 10000);
	m_paramsFile.addVariable(CFG_PRES_TIME_6, 10000);

	m_paramsFile.addVariable(CFG_RESP_TIME, 5000);
	m_paramsFile.addVariable(CFG_LOAD_TIME, 0);
	m_paramsFile.addVariable(CFG_FIX_TIME, 2000);
	m_paramsFile.addVariable(CFG_FIX_TIME_RNG, 2000);

	// Limiters
	m_paramsFile.addVariable(CFG_CONTRAST_LEV, "1");
	m_paramsFile.addVariable(CFG_PRES_TIME_LEV, "1");
	m_paramsFile.addVariable(CFG_FIX_LOC_NO, 8);
	m_paramsFile.addVariable(CFG_NUMB_TRIAL, 50);
	m_paramsFile.addVariable(CFG_NUMB_IMG, 50);

	// Exp Params
	//m_paramsFile.addVariable(CFG_SPAT_FREQ, 0);
	//m_paramsFile.addVariable(CFG_ECC, 0.f);
	m_paramsFile.addVariable(CFG_SPAT_FREQ_LEV, "2");
	m_paramsFile.addVariable(CFG_SPAT_FREQ_0, "1");
	m_paramsFile.addVariable(CFG_SPAT_FREQ_1, "1");
	m_paramsFile.addVariable(CFG_CONTRAST_BACK, "1");


	m_paramsFile.addVariable(CFG_FIX_ECC, 250);
	m_paramsFile.addVariable(CFG_FIX_AREA_C, 50);
	m_paramsFile.addVariable(CFG_FIX_AREA_P, 50);
	m_paramsFile.addVariable(CFG_FIX_ANGLE_OFFSET, 0);

	m_paramsFile.addVariable(CFG_CONTRAST_0, "0");
	m_paramsFile.addVariable(CFG_CONTRAST_1, "1");
	m_paramsFile.addVariable(CFG_CONTRAST_2, "2");
	m_paramsFile.addVariable(CFG_CONTRAST_3, "3");
	m_paramsFile.addVariable(CFG_CONTRAST_4, "4");
	m_paramsFile.addVariable(CFG_CONTRAST_5, "5");
	m_paramsFile.addVariable(CFG_CONTRAST_6, "6");

	m_paramsFile.addVariable(CFG_ECC_LEV, "2"); 
	m_paramsFile.addVariable(CFG_ECC_0, "0");
	m_paramsFile.addVariable(CFG_ECC_1, "1");
	m_paramsFile.addVariable(CFG_ECC_2, "2");
	m_paramsFile.addVariable(CFG_ECC_3, "3");
	m_paramsFile.addVariable(CFG_ECC_4, "4");
	m_paramsFile.addVariable(CFG_ECC_5, "5");
	m_paramsFile.addVariable(CFG_ECC_6, "6");

	m_paramsFile.addVariable(CFG_DATA_DESTINATION, std::string(""));

	m_paramsFile.addVariable(CFG_DEBUG, 1);
	m_paramsFile.addVariable(CFG_NRECAL, 1);
	m_paramsFile.addVariable(CFG_AWID, 1.f);
	m_paramsFile.addVariable(CFG_GSD, 1.f);
	m_paramsFile.addVariable(CFG_GSD2, 1.f);
	m_paramsFile.addVariable(CFG_RAMP_TIME, 5.0f);
	m_paramsFile.addVariable(CFG_GT, std::string(""));
}


// The one and only CEyeChartEApp object

CPostSaccSens_EccApp theApp;


// CEyeChartEApp initialization

BOOL CPostSaccSens_EccApp::InitInstance()
{
	INITCOMMONCONTROLSEX InitCtrls;
	InitCtrls.dwSize = sizeof(InitCtrls);
	InitCtrls.dwICC = ICC_WIN95_CLASSES;
	InitCommonControlsEx(&InitCtrls);
	CWinApp::InitInstance();
	AfxEnableControlContainer();
	CShellManager *pShellManager = new CShellManager;
	CMFCVisualManager::SetDefaultManager(RUNTIME_CLASS(CMFCVisualManagerWindows));
	SetRegistryKey(_T("EyeRIS @ APLab"));

	// EyeRIS ------------------------------------------------------------------
	// Initialize the library
	// Insert this line if the library needs to load a specific configuration file
	CWinEMIL::Instance()->initialize("emil-library-monocular-singleROG.cfg");
	//CWinEMIL::Instance()->initialize("C:/EyeRIS/system/config/emil-library.cfg");
	// Put your code here ======================================================

	// Load the specified configuration file
	m_paramsFile.loadFile("Params.cfg");
	
	std::string DestinationDir = m_paramsFile.getDirectory(CFG_DATA_DESTINATION) +
		m_paramsFile.getString(CFG_SUBJECT_NAME) + "/calibration";

	int HResolution = m_paramsFile.getInteger(CFG_SCREEN_W);
	int VResolution = m_paramsFile.getInteger(CFG_SCREEN_H);
	int RefreshRate = m_paramsFile.getInteger(CFG_SCREEN_R);

	
	// to set the initial target window
	ExTarget* target = new ExTarget(HResolution, VResolution, RefreshRate);
	target->setBackgroundColor(127, 127, 127);
	CWinEMIL::Instance()->addExperiment(target);

	
	// updated call for eyeris 2.4.1 version
	if (0)
	{
		myAutoCalibration* autoCalibrator = new myAutoCalibration(HResolution, VResolution, RefreshRate);
		autoCalibrator->setOutputDir(DestinationDir);
		autoCalibrator->setBackgroundColor(127, 127, 127);
		autoCalibrator->setPointOffset(380, 380); //Set calibration grid size here, the parameters define the spacing between the calibration points
		CWinEMIL::Instance()->addExperiment(autoCalibrator);


		//////// right eye calibration
		myManualCalibrator2* manualCalibrator2 = new myManualCalibrator2(HResolution, VResolution, RefreshRate, EOS_EYE_1, false);
		manualCalibrator2->setOutputDir(DestinationDir);
		manualCalibrator2->setBackgroundColor(127, 127, 127);
		manualCalibrator2->setPointOffset(380, 380); //Set calibration grid size here, the parameters define the spacing between the calibration points
		CWinEMIL::Instance()->addExperiment(manualCalibrator2);
	}

	CWinEMIL::Instance()->addExperiment(new ExperimentBody(HResolution, VResolution, RefreshRate, &m_paramsFile));

	// =========================================================================
	//CEyeChartEDlg dlg;
	CEMILConsole dlg("emil-console-monocular-singleROG.cfg");
	//CEMILConsole dlg("C:/EyeRIS/system/config/emil-console.cfg");
	m_pMainWnd = &dlg;
	INT_PTR nResponse = dlg.DoModal();
	 
	// Destroy the the library
	CWinEMIL::Destroy();
	// EyeRIS ------------------------------------------------------------------
	
	// Delete the shell manager created above.
	if (pShellManager != nullptr)
	{
		delete pShellManager;
	}

#if !defined(_AFXDLL) && !defined(_AFX_NO_MFC_CONTROLS_IN_DIALOGS)
	ControlBarCleanUp();
#endif

	// Since the dialog has been closed, return FALSE so that we exit the
	//  application, rather than start the application's message pump.
	return FALSE;
}

