echo part_1_sim1_50.R

Rscript.exe part_1_sim1_50.R
Write-Host "Congratulations! Part 1 of simulation 1 executed successfully. Starting with part 2."

echo part_ff_sim1_50.R

for ($i=2; $i -lt 11; $i++){
    	Rscript.exe part_ff_sim1_50.R

	if($i -lt 10){
	Write-Host "Congratulations! Part $($i) of simulation 1 executed successfully. Starting with part $($i+1)."
	}
		else {
	Write-Host "Congratulations! Part $($i) of simulation 1 executed successfully. The simulation ended successfully."
	}
	
}