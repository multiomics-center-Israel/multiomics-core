' ============================================================
' multiomics-core — Launch Config Wizard (hidden terminal)
' Double-click to open the wizard in your browser.
' The black terminal window is hidden — only the browser opens.
' ============================================================
Dim WshShell, fso, scriptDir
Set WshShell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
scriptDir = fso.GetParentFolderName(WScript.ScriptFullName)
WshShell.Run "cmd /c cd /d """ & scriptDir & """ && Rscript run.R --wizard", 0, False
