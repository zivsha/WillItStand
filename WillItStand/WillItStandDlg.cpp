
// WillItStandDlg.cpp : implementation file
//

#include "stdafx.h"
#include "WillItStand.h"
#include "WillItStandDlg.h"
#include "afxdialogex.h"
#include <string>

#include "vec3.h"
#include "MeshModel.h"
#include "WillItStandEngine.h"
#include "Qhull.h"
#include <fstream>

#include <cstdio>   /* for printf() of help message */
#include <ostream>
#include <conio.h>

using std::cerr;
using std::cin;
using std::cout;
using std::endl;
using std::vector;

using namespace orgQhull;
#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CWillItStandDlg dialog



CWillItStandDlg::CWillItStandDlg(CWnd* pParent /*=NULL*/)
	: CDialogEx(CWillItStandDlg::IDD, pParent)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
    m_exportCH = false;
    m_exportCOM = false;
    m_exportResults = false;
}

void CWillItStandDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
    DDX_Text(pDX, IDC_EDIT_FILE_NAME, m_inputFilePath);
    DDX_Text(pDX, IDC_STATIC_RESULTS, m_output);
    DDX_Check(pDX, IDC_CHECK_CH, m_exportCH);
    DDX_Check(pDX, IDC_CHECK_RESULT, m_exportResults);
    DDX_Check(pDX, IDC_CHECK_COM, m_exportCOM);
}

BEGIN_MESSAGE_MAP(CWillItStandDlg, CDialogEx)
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON_INPUT_FILE, &CWillItStandDlg::OnBnClickedButtonInputFile)
    ON_BN_CLICKED(IDC_CHECK_RESULT, &CWillItStandDlg::OnBnClickedCheckResult)
    ON_BN_CLICKED(IDC_CHECK_CH, &CWillItStandDlg::OnBnClickedCheckCh)
    ON_BN_CLICKED(IDC_BUTTON_RUN, &CWillItStandDlg::OnBnClickedButtonRun)
    ON_BN_CLICKED(IDC_BUTTON_SAVE_RESULT_OUTPUT, &CWillItStandDlg::OnBnClickedButtonSaveResultOutput)
    ON_BN_CLICKED(IDC_CHECK_COM, &CWillItStandDlg::OnBnClickedCheckCom)
    ON_BN_CLICKED(IDC_BUTTON_OPEN_DROP, &CWillItStandDlg::OnBnClickedButtonOpenDrop)
END_MESSAGE_MAP()


// CWillItStandDlg message handlers

BOOL CWillItStandDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon

	// TODO: Add extra initialization here

	return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CWillItStandDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CWillItStandDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CWillItStandDlg::OnBnClickedButtonInputFile()
{
    TCHAR szFilters[] = _T("OBJ  Files (*.obj)|*.obj|All Files (*.*)|*.*||");

    CFileDialog dlg(TRUE, L"obj", L"", OFN_FILEMUSTEXIST | OFN_HIDEREADONLY, szFilters);
    dlg.m_ofn.lpstrTitle = L"Load .obj input file";
    if (dlg.DoModal() == IDOK) {
        m_inputFilePath = dlg.GetPathName();		// Full path and filename
        m_inputFileTitle = dlg.GetFileTitle();
        m_inputFileDir = dlg.GetFolderPath();
        m_hasRun = false;
        UpdateData(FALSE);
    }
}


void CWillItStandDlg::OnBnClickedCheckResult()
{
    m_exportResults = !m_exportResults;
    UpdateData(0);
}


void CWillItStandDlg::OnBnClickedCheckCh()
{
    m_exportCH = !m_exportCH;
    UpdateData(0);
}

bool CWillItStandDlg::IsOutputCH()
{
    return m_exportCH ? true: false;
}
bool CWillItStandDlg::IsOutputResult(){
    return m_exportResults ? true : false;
}
bool CWillItStandDlg::IsOutputCOM()
{
    return m_exportCOM ? true : false;
}
std::string CWillItStandDlg::GetFileName()
{
    CT2CA pszConvertedAnsiString(m_inputFilePath);
    return std::string(pszConvertedAnsiString);
}

void CWillItStandDlg::OnBnClickedButtonRun()
{
    if (m_inputFilePath.IsEmpty() || m_hasRun)
    {
        return;
    }
    m_hasRun = true;
    CT2CA pszConvertedAnsiString(m_inputFilePath);
    string fileName = string(pszConvertedAnsiString);
    CT2CA directory(m_inputFileDir);
    string dir = string(directory);
    WillItStand proj(fileName, m_exportResults ? true : false, m_exportCH ? true : false, m_exportCOM ? true : false, dir);
    ostringstream os;
    HCURSOR hWaitCursor = LoadCursor(NULL, IDC_APPSTARTING);
    if (hWaitCursor)
    {
        SetCursor(hWaitCursor);
    }
    bool result = proj.Run(os);

    HCURSOR hArrowCursor = LoadCursor(NULL, IDC_ARROW);
    if (hArrowCursor)
    {
        SetCursor(hArrowCursor);
    }
    theApp.m_pMainWnd->MessageBox(result ? CString("Finished processing model: ") + m_inputFileTitle : L"Error while running, see output for more info",
                                  result ? L"Process Complete" : L"Error", 
                                  result ? MB_ICONINFORMATION : MB_ICONWARNING);

    m_output = os.str().c_str();
    UpdateData(0);
}
void CWillItStandDlg::OnBnClickedButtonSaveResultOutput()
{
    if (m_inputFilePath.IsEmpty() || !m_hasRun || m_output.IsEmpty())
    {
        return;
    }
    TCHAR szFilters[] = _T("Text Files (*.txt)|*.txt|All Files (*.*)|*.*||");
    CFileDialog dlg(FALSE, L".txt", m_inputFileTitle + "_Results", OFN_HIDEREADONLY | OFN_OVERWRITEPROMPT, szFilters);
    dlg.m_ofn.lpstrTitle = L"Save Results to file";
    if (dlg.DoModal() == IDOK)
    {
        CT2CA caFileName(dlg.GetPathName());
        string fileName = string(caFileName);
        std::ofstream outFile;
        outFile.open(fileName.c_str());
        CT2CA caOutput(m_output);
        outFile << string(caOutput);
        outFile.close();
    }
}
void CWillItStandDlg::OnBnClickedCheckCom()
{
    m_exportCOM = !m_exportCOM;
}
void CWillItStandDlg::OnBnClickedButtonOpenDrop()
{
   if (m_hasRun && (m_exportCH || m_exportCOM || m_exportResults))
   {
       ShellExecute(NULL, L"open", m_inputFileDir, NULL, NULL, SW_SHOWDEFAULT);
   }
}
