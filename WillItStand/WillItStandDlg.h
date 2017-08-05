
// WillItStandDlg.h : header file
//

#pragma once

#include <string>

// CWillItStandDlg dialog
class CWillItStandDlg : public CDialogEx
{
// Construction
    CString m_inputFilePath;
    CString m_inputFileTitle;
    CString m_inputFileDir;
    CString m_output;
    BOOL m_exportCH;
    BOOL m_exportCOM;
    BOOL m_exportResults;
    BOOL m_hasRun;
public:
	CWillItStandDlg(CWnd* pParent = NULL);	// standard constructor
    bool IsOutputCH();
    bool IsOutputResult();
    bool IsOutputCOM();
    std::string GetFileName(); 
    
// Dialog Data
	enum { IDD = IDD_WILLITSTAND_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
    afx_msg void OnBnClickedButtonInputFile();
    afx_msg void OnBnClickedCheckResult();
    afx_msg void OnBnClickedCheckCh();
    afx_msg void OnBnClickedCheckIntersecting();
    afx_msg void OnBnClickedButtonRun();
    afx_msg void OnBnClickedButtonSaveResultOutput();
    afx_msg void OnBnClickedCheckCom();
    afx_msg void OnBnClickedButtonOpenDrop();
};
