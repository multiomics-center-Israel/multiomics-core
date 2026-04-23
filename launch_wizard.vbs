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

' Run the wizard with hidden terminal (0 = hidden, False = don't wait)
WshShell.Run "cmd /c cd /d """ & scriptDir & """ && " & setPandoc & "Rscript --vanilla run.R --wizard", 0, False
