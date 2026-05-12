' Creates a desktop shortcut for the multiomics-core wizard
Dim WshShell, fso, scriptDir, shortcut, desktopPath
Set WshShell = CreateObject("WScript.Shell")
Set fso = CreateObject("Scripting.FileSystemObject")
scriptDir = fso.GetParentFolderName(WScript.ScriptFullName)
desktopPath = WshShell.SpecialFolders("Desktop")

Set shortcut = WshShell.CreateShortcut(desktopPath & "\Multiomics Pipeline.lnk")
shortcut.TargetPath = scriptDir & "\launch_wizard.vbs"
shortcut.WorkingDirectory = scriptDir
shortcut.IconLocation = "shell32.dll,13"
shortcut.Description = "Launch multiomics-core Config Wizard"
shortcut.Save

WScript.Echo "[OK] Desktop shortcut created: Multiomics Pipeline"
