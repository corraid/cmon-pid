$kst = "..\Kst2\bin\kst2.exe"
$imgs = "images\"
New-Item -Path $imgs -Type Directory -Force
$csvfiles = Get-ChildItem "." -Filter *.csv
foreach ($f in $csvfiles) {
	$png = $imgs + $f.Basename + ".png"
	&$kst sim_fast_sample.kst -F $f --pngHeight 1500 --pngWidth 1500 --png $png
}