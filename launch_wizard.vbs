' ============================================================
' multiomics-core - Launch Config Wizard (hidden terminal)
' Double-click to open the wizard in your browser.
' ============================================================
Dim WshShell, fso, scriptDir
Set WshShell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
scriptDir = fso.GetParentFolderName(WScript.ScriptFullName)

' Show a brief notification so user knows it's loading
WshShell.Popup "Starting multiomics-core wizard..." & vbCrLf & vbCrLf & "Your browser will open shortly. This window will close automatically.", 30, "Multiomics Pipeline", 64

' Run the wizard with hidden terminal (0 = hidden, False = don't wait)
WshShell.Run "cmd /c cd /d """ & scriptDir & """ && Rscript run.R --wizard", 0, False
