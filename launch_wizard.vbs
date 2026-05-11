' ============================================================
' multiomics-core - Launch Config Wizard (hidden terminal)
' Double-click to open the wizard in your browser.
' ============================================================
Dim WshShell, fso, scriptDir
Set WshShell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
scriptDir = fso.GetParentFolderName(WScript.ScriptFullName)

' Show a brief notification so user knows it's loading
WshShell.Popup "Starting multiomics-core wizard..." & vbCrLf & vbCrLf & "Your browser will open shortly." & vbCrLf & "If it doesn't open within ~30s, open wizard_launch.log in this folder for details.", 20, "Multiomics Pipeline", 64

' If the bundled portable pandoc is present, point rmarkdown at it via
' RSTUDIO_PANDOC. This avoids pandoc "error 127" when the install path
' contains spaces or when RStudio is not installed on the machine.
Dim pandocDir
pandocDir = scriptDir & "\tools\pandoc"
Dim setPandoc
setPandoc = ""
If fso.FileExists(pandocDir & "\pandoc.exe") Then
    setPandoc = "set ""RSTUDIO_PANDOC=" & pandocDir & """ && "
End If

' Run the wizard with hidden terminal (0 = hidden, False = don't wait).
' Tee output into wizard_launch.log so failures aren't silent — if the
' browser never opens, the user can open that file to see the error.
Dim logPath
logPath = scriptDir & "\wizard_launch.log"
WshShell.Run "cmd /c cd /d """ & scriptDir & """ && " & setPandoc & "Rscript --vanilla run.R --wizard > """ & logPath & """ 2>&1", 0, False
