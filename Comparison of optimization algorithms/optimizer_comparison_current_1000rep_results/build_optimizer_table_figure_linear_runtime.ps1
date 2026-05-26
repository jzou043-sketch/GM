Add-Type -AssemblyName System.Drawing
Add-Type -AssemblyName System.IO.Compression.FileSystem

$ErrorActionPreference = "Stop"

$base = "C:\Users\15280\OneDrive\Desktop\optimizer_comparison_current_1000rep\optimizer_comparison_current_1000rep_results"
$summaryPath = Join-Path $base "benchmark_1000rep_optimizer_current_summary.csv"
$rawPath = Join-Path $base "benchmark_1000rep_optimizer_current_replicate_best.csv"
$outDir = "C:\Users\15280\OneDrive\Desktop\optimizer_algorithm_comparison_log_runtime_outputs"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$order = @(
  "Particle swarm optimization",
  "Tabu search",
  "Proposed method",
  "Hill climbing",
  "Simulated annealing",
  "Genetic algorithm",
  "Differential evolution"
)

$summary = Import-Csv -LiteralPath $summaryPath
$raw = Import-Csv -LiteralPath $rawPath

function ToD($x) {
  return [double]::Parse([string]$x, [System.Globalization.CultureInfo]::InvariantCulture)
}

function F3($x) {
  return ([double]$x).ToString("0.000", [System.Globalization.CultureInfo]::InvariantCulture)
}

function F1($x) {
  return ([double]$x).ToString("0.0", [System.Globalization.CultureInfo]::InvariantCulture)
}

function XmlEscape([string]$s) {
  return [System.Security.SecurityElement]::Escape($s)
}

# ---------------------------------------------------------
# Three-line table data
# ---------------------------------------------------------
$tableRows = foreach ($alg in $order) {
  $r = $summary | Where-Object { $_.Algorithm -eq $alg } | Select-Object -First 1
  [pscustomobject]@{
    Algorithm = $alg
    "Mean best objective" = F3 (ToD $r.MeanBestObjective)
    SD = F3 (ToD $r.SDBestObjective)
    "Mean runtime (s)" = F3 (ToD $r.MeanRuntimeSec)
    "Top-hit rate (%)" = F1 (ToD $r.TopHitRatePct)
    "Mean regret" = F3 (ToD $r.MeanRegret)
  }
}

$tableCsv = Join-Path $outDir "optimizer_algorithm_summary_three_line_table.csv"
$tableRows | Export-Csv -LiteralPath $tableCsv -NoTypeInformation -Encoding UTF8

# ---------------------------------------------------------
# Three-line table DOCX
# ---------------------------------------------------------
function CellXml([string]$text, [int]$width, [bool]$bold = $false, [string]$jc = "center", [bool]$bottom = $false) {
  $b = if ($bold) { "<w:b/>" } else { "" }
  $border = if ($bottom) {
    '<w:tcBorders><w:bottom w:val="single" w:sz="8" w:space="0" w:color="000000"/></w:tcBorders>'
  } else {
    ""
  }
  $escaped = XmlEscape $text
  return @"
<w:tc>
  <w:tcPr><w:tcW w:w="$width" w:type="dxa"/>$border<w:tcMar><w:top w:w="90" w:type="dxa"/><w:left w:w="90" w:type="dxa"/><w:bottom w:w="90" w:type="dxa"/><w:right w:w="90" w:type="dxa"/></w:tcMar></w:tcPr>
  <w:p><w:pPr><w:jc w:val="$jc"/></w:pPr><w:r><w:rPr>$b<w:sz w:val="20"/><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/></w:rPr><w:t>$escaped</w:t></w:r></w:p>
</w:tc>
"@
}

$headers = @("Algorithm", "Mean best objective", "SD", "Mean runtime (s)", "Top-hit rate (%)", "Mean regret")
$widths = @(3600, 1900, 1100, 1700, 1600, 1400)
$rowsXml = New-Object System.Text.StringBuilder

[void]$rowsXml.Append("<w:tr>")
for ($i = 0; $i -lt $headers.Count; $i++) {
  [void]$rowsXml.Append((CellXml $headers[$i] $widths[$i] $true "center" $true))
}
[void]$rowsXml.Append("</w:tr>")

foreach ($r in $tableRows) {
  [void]$rowsXml.Append("<w:tr>")
  [void]$rowsXml.Append((CellXml $r.Algorithm $widths[0] $false "left" $false))
  [void]$rowsXml.Append((CellXml $r."Mean best objective" $widths[1]))
  [void]$rowsXml.Append((CellXml $r.SD $widths[2]))
  [void]$rowsXml.Append((CellXml $r."Mean runtime (s)" $widths[3]))
  [void]$rowsXml.Append((CellXml $r."Top-hit rate (%)" $widths[4]))
  [void]$rowsXml.Append((CellXml $r."Mean regret" $widths[5]))
  [void]$rowsXml.Append("</w:tr>")
}

$gridXml = ($widths | ForEach-Object { '<w:gridCol w:w="' + $_ + '"/>' }) -join "`n"
$tableBodyXml = $rowsXml.ToString()

$documentXml = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<w:document xmlns:w="http://schemas.openxmlformats.org/wordprocessingml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
  <w:body>
    <w:p><w:pPr><w:jc w:val="center"/></w:pPr><w:r><w:rPr><w:b/><w:sz w:val="24"/><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/></w:rPr><w:t>Table. Optimization-algorithm comparison in the simulated genomic mating problem</w:t></w:r></w:p>
    <w:tbl>
      <w:tblPr>
        <w:tblW w:w="0" w:type="auto"/>
        <w:jc w:val="center"/>
        <w:tblBorders>
          <w:top w:val="single" w:sz="12" w:space="0" w:color="000000"/>
          <w:left w:val="nil"/>
          <w:bottom w:val="single" w:sz="12" w:space="0" w:color="000000"/>
          <w:right w:val="nil"/>
          <w:insideH w:val="nil"/>
          <w:insideV w:val="nil"/>
        </w:tblBorders>
        <w:tblLook w:firstRow="1" w:lastRow="0" w:firstColumn="0" w:lastColumn="0" w:noHBand="1" w:noVBand="1"/>
      </w:tblPr>
      <w:tblGrid>
        $gridXml
      </w:tblGrid>
      $tableBodyXml
    </w:tbl>
    <w:p><w:r><w:rPr><w:sz w:val="18"/><w:rFonts w:ascii="Times New Roman" w:hAnsi="Times New Roman"/></w:rPr><w:t>Note. Top-hit rate was calculated within each simulation replicate relative to the best objective observed among all algorithms. Regret was calculated as the replicate-specific best observed objective minus the objective obtained by each algorithm.</w:t></w:r></w:p>
    <w:sectPr><w:pgSz w:w="16838" w:h="11906" w:orient="landscape"/><w:pgMar w:top="900" w:right="900" w:bottom="900" w:left="900" w:header="720" w:footer="720" w:gutter="0"/></w:sectPr>
  </w:body>
</w:document>
"@

$contentTypes = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
  <Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
  <Default Extension="xml" ContentType="application/xml"/>
  <Override PartName="/word/document.xml" ContentType="application/vnd.openxmlformats-officedocument.wordprocessingml.document.main+xml"/>
  <Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
  <Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>
"@

$rels = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
  <Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="word/document.xml"/>
</Relationships>
"@

$docRels = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships"/>
"@

$core = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:title>Optimizer comparison table</dc:title><dc:creator>Codex</dc:creator></cp:coreProperties>
"@

$app = @"
<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"><Application>Microsoft Word</Application></Properties>
"@

$tmpDocxDir = Join-Path $env:TEMP ("docx_" + [Guid]::NewGuid().ToString("N"))
New-Item -ItemType Directory -Path $tmpDocxDir | Out-Null
New-Item -ItemType Directory -Path (Join-Path $tmpDocxDir "_rels") | Out-Null
New-Item -ItemType Directory -Path (Join-Path $tmpDocxDir "word") | Out-Null
New-Item -ItemType Directory -Path (Join-Path $tmpDocxDir "word\_rels") | Out-Null
New-Item -ItemType Directory -Path (Join-Path $tmpDocxDir "docProps") | Out-Null

Set-Content -LiteralPath (Join-Path $tmpDocxDir "[Content_Types].xml") -Value $contentTypes -Encoding UTF8
Set-Content -LiteralPath (Join-Path $tmpDocxDir "_rels\.rels") -Value $rels -Encoding UTF8
Set-Content -LiteralPath (Join-Path $tmpDocxDir "word\document.xml") -Value $documentXml -Encoding UTF8
Set-Content -LiteralPath (Join-Path $tmpDocxDir "word\_rels\document.xml.rels") -Value $docRels -Encoding UTF8
Set-Content -LiteralPath (Join-Path $tmpDocxDir "docProps\core.xml") -Value $core -Encoding UTF8
Set-Content -LiteralPath (Join-Path $tmpDocxDir "docProps\app.xml") -Value $app -Encoding UTF8

$tableDocx = Join-Path $outDir "optimizer_algorithm_summary_three_line_table.docx"
if (Test-Path $tableDocx) { Remove-Item -LiteralPath $tableDocx -Force }
[System.IO.Compression.ZipFile]::CreateFromDirectory($tmpDocxDir, $tableDocx)
Remove-Item -LiteralPath $tmpDocxDir -Recurse -Force

# ---------------------------------------------------------
# Boxplot figure, runtime uses log10 scale
# ---------------------------------------------------------
function Percentile([double[]]$arr, [double]$p) {
  $a = @($arr | Sort-Object)
  $n = $a.Count
  if ($n -eq 0) { return [double]::NaN }
  if ($n -eq 1) { return [double]$a[0] }
  $pos = ($n - 1) * $p
  $lo = [math]::Floor($pos)
  $hi = [math]::Ceiling($pos)
  if ($lo -eq $hi) { return [double]$a[$lo] }
  return [double]$a[$lo] + ($pos - $lo) * ([double]$a[$hi] - [double]$a[$lo])
}

function BoxStats([double[]]$x) {
  $q1 = Percentile $x 0.25
  $med = Percentile $x 0.50
  $q3 = Percentile $x 0.75
  $iqr = $q3 - $q1
  $loFence = $q1 - 1.5 * $iqr
  $hiFence = $q3 + 1.5 * $iqr
  $in = @($x | Where-Object { $_ -ge $loFence -and $_ -le $hiFence })
  [pscustomobject]@{
    Q1 = $q1
    Median = $med
    Q3 = $q3
    Low = ($in | Measure-Object -Minimum).Minimum
    High = ($in | Measure-Object -Maximum).Maximum
    Outliers = @($x | Where-Object { $_ -lt $loFence -or $_ -gt $hiFence })
  }
}

function NiceTicks([double]$min, [double]$max, [int]$n = 6) {
  if ($max -le $min) { return @($min, $max) }
  $range = $max - $min
  $rawStep = $range / [math]::Max(1, ($n - 1))
  $pow = [math]::Pow(10, [math]::Floor([math]::Log10($rawStep)))
  $steps = @(1, 2, 5, 10)
  $step = $pow
  foreach ($s in $steps) {
    if ($rawStep -le $s * $pow) {
      $step = $s * $pow
      break
    }
  }
  $start = [math]::Floor($min / $step) * $step
  $end = [math]::Ceiling($max / $step) * $step
  $ticks = @()
  $v = $start
  while ($v -le $end + 1e-9) {
    $ticks += [double]$v
    $v += $step
  }
  return $ticks
}

function DrawRotatedText($g, [string]$text, $font, $brush, [float]$x, [float]$y, [float]$angle) {
  $state = $g.Save()
  $g.TranslateTransform($x, $y)
  $g.RotateTransform($angle)
  $g.DrawString($text, $font, $brush, 0, 0)
  $g.Restore($state)
}

$plotData = @{}
foreach ($alg in $order) {
  $rows = @($raw | Where-Object { $_.Algorithm -eq $alg })
  $plotData[$alg] = [pscustomobject]@{
    Obj = [double[]]($rows | ForEach-Object { ToD $_.BestObjective })
    Rt = [double[]]($rows | ForEach-Object { ToD $_.TotalElapsedSec })
  }
}

$width = 9000
$height = 4300
$bmp = New-Object System.Drawing.Bitmap $width, $height
$bmp.SetResolution(1000, 1000)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::AntiAliasGridFit
$g.Clear([System.Drawing.Color]::White)

$fontTitle = New-Object System.Drawing.Font("Arial", 92, [System.Drawing.FontStyle]::Regular, [System.Drawing.GraphicsUnit]::Pixel)
$fontAxis = New-Object System.Drawing.Font("Arial", 70, [System.Drawing.FontStyle]::Regular, [System.Drawing.GraphicsUnit]::Pixel)
$fontTick = New-Object System.Drawing.Font("Arial", 56, [System.Drawing.FontStyle]::Regular, [System.Drawing.GraphicsUnit]::Pixel)
$fontLabel = New-Object System.Drawing.Font("Arial", 52, [System.Drawing.FontStyle]::Regular, [System.Drawing.GraphicsUnit]::Pixel)
$black = [System.Drawing.Brushes]::Black
$grayPen = New-Object System.Drawing.Pen([System.Drawing.Color]::FromArgb(225,225,225), 5)
$axisPen = New-Object System.Drawing.Pen([System.Drawing.Color]::Black, 6)
$thinPen = New-Object System.Drawing.Pen([System.Drawing.Color]::Black, 5)
$outBrush = New-Object System.Drawing.SolidBrush([System.Drawing.Color]::FromArgb(55,55,55))

$colors = @(
  [System.Drawing.Color]::FromArgb(75,58,115),
  [System.Drawing.Color]::FromArgb(68,87,130),
  [System.Drawing.Color]::FromArgb(46,121,131),
  [System.Drawing.Color]::FromArgb(43,139,132),
  [System.Drawing.Color]::FromArgb(77,171,115),
  [System.Drawing.Color]::FromArgb(102,180,94),
  [System.Drawing.Color]::FromArgb(167,190,68)
)

function DrawPanel($g, [System.Drawing.RectangleF]$rect, [string]$title, [string]$ylabel, [string]$type) {
  $titleSize = $g.MeasureString($title, $fontTitle)
  $g.DrawString($title, $fontTitle, $black, $rect.X + $rect.Width / 2 - $titleSize.Width / 2, $rect.Y - 155)

  $left = $rect.X + 430
  $right = $rect.X + $rect.Width - 130
  $top = $rect.Y + 80
  $bottom = $rect.Y + $rect.Height - 780
  $plotW = $right - $left
  $plotH = $bottom - $top

  if ($type -eq "obj") {
    $all = [double[]]($order | ForEach-Object { $plotData[$_].Obj })
    $dataMin = ($all | Measure-Object -Minimum).Minimum
    $dataMax = ($all | Measure-Object -Maximum).Maximum
    $ticks = NiceTicks ($dataMin - 8) ($dataMax + 8) 6
  } else {
    $allRaw = [double[]]($order | ForEach-Object { $plotData[$_].Rt | Where-Object { $_ -gt 0 } })
    $dataMin = ($allRaw | Measure-Object -Minimum).Minimum
    $dataMax = ($allRaw | Measure-Object -Maximum).Maximum
    $pMin = [math]::Floor([math]::Log10($dataMin))
    $pMax = [math]::Ceiling([math]::Log10($dataMax))
    $ticks = @()
    for ($p = $pMin; $p -le $pMax; $p++) {
      foreach ($m in @(1, 2, 5)) {
        $t = $m * [math]::Pow(10, $p)
        if ($t -ge $dataMin * 0.8 -and $t -le $dataMax * 1.2) {
          $ticks += [double]$t
        }
      }
    }
    if ($ticks.Count -lt 2) {
      $ticks = @([math]::Pow(10, $pMin), [math]::Pow(10, $pMax))
    }
  }

  $yMin = if ($type -eq "obj") { [double]($ticks | Select-Object -First 1) } else { [math]::Log10([double]($ticks | Select-Object -First 1)) }
  $yMax = if ($type -eq "obj") { [double]($ticks | Select-Object -Last 1) } else { [math]::Log10([double]($ticks | Select-Object -Last 1)) }

  foreach ($t in $ticks) {
    $tv = if ($type -eq "obj") { [double]$t } else { [math]::Log10([double]$t) }
    $yy = $bottom - (($tv - $yMin) / ($yMax - $yMin)) * $plotH
    $g.DrawLine($grayPen, [float]$left, [float]$yy, [float]$right, [float]$yy)
    $lab = if ($type -eq "obj") { ([double]$t).ToString("0") } else { ([double]$t).ToString("0") }
    $sz = $g.MeasureString($lab, $fontTick)
    $g.DrawString($lab, $fontTick, $black, $left - $sz.Width - 30, [float]$yy - $sz.Height / 2)
  }

  $g.DrawLine($axisPen, [float]$left, [float]$top, [float]$left, [float]$bottom)
  $g.DrawLine($axisPen, [float]$left, [float]$bottom, [float]$right, [float]$bottom)
  DrawRotatedText $g $ylabel $fontAxis $black ($rect.X + 30) ($top + $plotH / 2 + 260) -90

  $n = $order.Count
  $step = $plotW / $n
  $boxW = [math]::Min(330, $step * 0.62)

  for ($i = 0; $i -lt $n; $i++) {
    $alg = $order[$i]
    $x = $left + $step * ($i + 0.5)
    $vals = if ($type -eq "obj") {
      $plotData[$alg].Obj
    } else {
      [double[]]($plotData[$alg].Rt | Where-Object { $_ -gt 0 } | ForEach-Object { [math]::Log10([double]$_) })
    }
    $bs = BoxStats $vals

    $map = {
      param($v)
      $bottom - (($v - $yMin) / ($yMax - $yMin)) * $plotH
    }

    $yQ1 = & $map $bs.Q1
    $yQ3 = & $map $bs.Q3
    $yMed = & $map $bs.Median
    $yLow = & $map $bs.Low
    $yHigh = & $map $bs.High

    $fill = New-Object System.Drawing.SolidBrush($colors[$i])
    $g.DrawLine($thinPen, [float]$x, [float]$yHigh, [float]$x, [float]$yQ3)
    $g.DrawLine($thinPen, [float]$x, [float]$yQ1, [float]$x, [float]$yLow)
    $g.DrawLine($thinPen, [float]($x - $boxW * 0.28), [float]$yHigh, [float]($x + $boxW * 0.28), [float]$yHigh)
    $g.DrawLine($thinPen, [float]($x - $boxW * 0.28), [float]$yLow, [float]($x + $boxW * 0.28), [float]$yLow)
    $rectBox = New-Object System.Drawing.RectangleF ([float]($x - $boxW / 2)), ([float][math]::Min($yQ1, $yQ3)), ([float]$boxW), ([float][math]::Abs($yQ3 - $yQ1))
    $g.FillRectangle($fill, $rectBox)
    $g.DrawRectangle($thinPen, $rectBox.X, $rectBox.Y, $rectBox.Width, $rectBox.Height)
    $g.DrawLine($axisPen, [float]($x - $boxW / 2), [float]$yMed, [float]($x + $boxW / 2), [float]$yMed)

    $outs = @($bs.Outliers)
    if ($outs.Count -gt 80) {
      $loOut = @($outs | Sort-Object | Select-Object -First 40)
      $hiOut = @($outs | Sort-Object | Select-Object -Last 40)
      $outs = @($loOut + $hiOut)
    }
    foreach ($o in $outs) {
      if ($o -ge $yMin -and $o -le $yMax) {
        $yy = & $map ([double]$o)
        $g.FillEllipse($outBrush, [float]($x - 8), [float]($yy - 8), 16, 16)
      }
    }

    $g.DrawLine($axisPen, [float]$x, [float]$bottom, [float]$x, [float]($bottom + 20))
    DrawRotatedText $g $alg $fontLabel $black ([float]($x - 135)) ([float]($bottom + 120)) -45
    $fill.Dispose()
  }
}

$panel1 = New-Object System.Drawing.RectangleF 60, 380, 4320, 3600
$panel2 = New-Object System.Drawing.RectangleF 4560, 380, 4320, 3600
DrawPanel $g $panel1 "A. Solution Quality (Objective Value)" "Objective Value" "obj"
DrawPanel $g $panel2 "B. Runtime (Seconds, Log Scale)" "Runtime (seconds)" "rt"

$figTiff = Join-Path $outDir "optimizer_algorithm_comparison_boxplots_log_runtime_1000dpi.tiff"
$figPng = Join-Path $outDir "optimizer_algorithm_comparison_boxplots_log_runtime_preview.png"
if (Test-Path $figTiff) { Remove-Item -LiteralPath $figTiff -Force }
if (Test-Path $figPng) { Remove-Item -LiteralPath $figPng -Force }

$codec = [System.Drawing.Imaging.ImageCodecInfo]::GetImageEncoders() | Where-Object { $_.MimeType -eq "image/tiff" } | Select-Object -First 1
$params = New-Object System.Drawing.Imaging.EncoderParameters 1
$params.Param[0] = New-Object System.Drawing.Imaging.EncoderParameter ([System.Drawing.Imaging.Encoder]::Compression), ([long][System.Drawing.Imaging.EncoderValue]::CompressionLZW)
$bmp.Save($figTiff, $codec, $params)
$bmp.Save($figPng, [System.Drawing.Imaging.ImageFormat]::Png)

$g.Dispose()
$bmp.Dispose()

Get-Item -LiteralPath $figTiff, $figPng, $tableCsv, $tableDocx | Select-Object FullName,Length,LastWriteTime
