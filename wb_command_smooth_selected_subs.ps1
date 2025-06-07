# Specify the root path to search
$rootPath = ""  # Change to your desired path

# Read the list of subnames to include from a file
$includeFilePath = "results\sub_ids_include.txt" 
if (-Not (Test-Path $includeFilePath)) {
    Write-Error "Include file not found: $includeFilePath"
    return
}

# Assume each line contains only a number (subname), trim whitespace, and store in an array
$includeSubs = Get-Content $includeFilePath | ForEach-Object { $_.Trim() } | Where-Object { $_ -ne "" }

# Output the loaded includeSubs list for confirmation
Write-Output "=== Start of includeSubs list (Total: $($includeSubs.Count)) ==="
$includeSubs | ForEach-Object { Write-Output $_ }
Write-Output "=== End of includeSubs list ==="
Write-Output ""

# Find folders directly under rootPath that contain "rfMRI" in their names
$data_folders = Get-ChildItem -Path $rootPath -Directory | Where-Object { $_.Name -match "rfMRI" }

foreach ($folder in $data_folders) {
    # Split folder name by '_'
    $parts = $folder.Name -split '_'
    
    # Expect at least 4 elements like in MATLAB code
    if ($parts.Length -ge 4) {
        # Extract the first and fourth elements (PowerShell uses 0-based index)
        $subname = $parts[0]
        $REST_num = $parts[3]

        # Check if subname is in the include list
        if (-Not ($includeSubs -contains $subname)) {
            Write-Output "Skipping subject '$subname' (not in the include list)"
            continue
        }

        # Check if REST_num is "REST1"
        if ($REST_num -ne "REST1") {
            Write-Output "Skipping subject '$subname' (REST_num is not 'REST1', but '$REST_num')"
            continue
        }
        
        # Construct LR file path
        $datafile_path_LR = Join-Path $folder.FullName $subname
        $datafile_path_LR = Join-Path $datafile_path_LR "MNINonLinear"
        $datafile_path_LR = Join-Path $datafile_path_LR "Results"
        $datafile_path_LR = Join-Path $datafile_path_LR ("rfMRI_" + $REST_num + "_LR")
        $input_datafile_path_LR = Join-Path $datafile_path_LR ("rfMRI_" + $REST_num + "_LR_Atlas_hp2000_clean.dtseries.nii")
        $output_datafile_path_LR = Join-Path $datafile_path_LR ("s6_rfMRI_" + $REST_num + "_LR_Atlas_hp2000_clean.dtseries.nii")
        
        # Construct RL file path
        $datafile_path_RL = Join-Path $folder.FullName $subname
        $datafile_path_RL = Join-Path $datafile_path_RL "MNINonLinear"
        $datafile_path_RL = Join-Path $datafile_path_RL "Results"
        $datafile_path_RL = Join-Path $datafile_path_RL ("rfMRI_" + $REST_num + "_RL")
        $input_datafile_path_RL = Join-Path $datafile_path_RL ("rfMRI_" + $REST_num + "_RL_Atlas_hp2000_clean.dtseries.nii")
        $output_datafile_path_RL = Join-Path $datafile_path_RL ("s6_rfMRI_" + $REST_num + "_RL_Atlas_hp2000_clean.dtseries.nii")

        # Check if smoothed RL output file already exists
        if (Test-Path $output_datafile_path_RL) {
            Write-Output "Skipping subject '$subname' (smoothing already completed)"
            continue
        }

        # Generate wb_command strings (wrap paths in quotes to handle spaces)
        $cmd_LR = "wb_command -cifti-smoothing `"$input_datafile_path_LR`" 6 6 COLUMN `"$output_datafile_path_LR`" -fwhm -left-surface C:\Users\YJ-second\Documents\ProjectsYJ\2025_DMD_with_ICA_fMRI\atlas\Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii -right-surface C:\Users\YJ-second\Documents\ProjectsYJ\2025_DMD_with_ICA_fMRI\atlas\Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii"
        $cmd_RL = "wb_command -cifti-smoothing `"$input_datafile_path_RL`" 6 6 COLUMN `"$output_datafile_path_RL`" -fwhm -left-surface C:\Users\YJ-second\Documents\ProjectsYJ\2025_DMD_with_ICA_fMRI\atlas\Q1-Q6_R440.L.midthickness.32k_fs_LR.surf.gii -right-surface C:\Users\YJ-second\Documents\ProjectsYJ\2025_DMD_with_ICA_fMRI\atlas\Q1-Q6_R440.R.midthickness.32k_fs_LR.surf.gii"
        
        # Output the commands for verification
        Write-Output "Subject: $subname, REST: $REST_num"
        Write-Output $cmd_LR
        Write-Output $cmd_RL
        Write-Output "------------------------------------"

        # Execute wb_command
        Invoke-Expression $cmd_LR
        Invoke-Expression $cmd_RL
    }
    else {
        Write-Output "Skipping folder '$($folder.FullName)': Folder name does not match expected format (at least 4 parts separated by '_')"
    }
}
