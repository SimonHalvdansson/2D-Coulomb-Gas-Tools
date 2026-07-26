param(
  [Parameter(Mandatory = $true)]
  [string]$OutputDirectory
)

$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Drawing

$maskSize = 384
$targetExtent = [single]($maskSize * 0.84)
$characters = @()
foreach ($code in (65..90 + 97..122 + 48..57)) {
  $characters += [char]$code
}

[System.IO.Directory]::CreateDirectory($OutputDirectory) | Out-Null
$fontFamily = [System.Drawing.FontFamily]::new('Arial Black')
$format = [System.Drawing.StringFormat]::GenericTypographic

foreach ($character in $characters) {
  $bitmap = [System.Drawing.Bitmap]::new(
    $maskSize,
    $maskSize,
    [System.Drawing.Imaging.PixelFormat]::Format32bppArgb
  )
  $graphics = [System.Drawing.Graphics]::FromImage($bitmap)
  $graphics.Clear([System.Drawing.Color]::Transparent)
  $graphics.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
  $graphics.PixelOffsetMode = [System.Drawing.Drawing2D.PixelOffsetMode]::HighQuality

  $path = [System.Drawing.Drawing2D.GraphicsPath]::new()
  $path.AddString(
    [string]$character,
    $fontFamily,
    [int][System.Drawing.FontStyle]::Regular,
    [single]344,
    [System.Drawing.PointF]::new(0, 0),
    $format
  )

  $bounds = $path.GetBounds()
  $scale = [Math]::Min($targetExtent / $bounds.Width, $targetExtent / $bounds.Height)
  $normalize = [System.Drawing.Drawing2D.Matrix]::new()
  $normalize.Translate(-$bounds.X, -$bounds.Y)
  $path.Transform($normalize)
  $normalize.Dispose()

  $scaleMatrix = [System.Drawing.Drawing2D.Matrix]::new()
  $scaleMatrix.Scale([single]$scale, [single]$scale)
  $path.Transform($scaleMatrix)
  $scaleMatrix.Dispose()

  $scaledBounds = $path.GetBounds()
  $center = [System.Drawing.Drawing2D.Matrix]::new()
  $center.Translate(
    [single](($maskSize - $scaledBounds.Width) / 2 - $scaledBounds.X),
    [single](($maskSize - $scaledBounds.Height) / 2 - $scaledBounds.Y)
  )
  $path.Transform($center)
  $center.Dispose()

  $graphics.FillPath([System.Drawing.Brushes]::White, $path)
  $rect = [System.Drawing.Rectangle]::new(0, 0, $maskSize, $maskSize)
  $bitmapData = $bitmap.LockBits(
    $rect,
    [System.Drawing.Imaging.ImageLockMode]::ReadOnly,
    [System.Drawing.Imaging.PixelFormat]::Format32bppArgb
  )
  $pixelBytes = [byte[]]::new([Math]::Abs($bitmapData.Stride) * $maskSize)
  [System.Runtime.InteropServices.Marshal]::Copy(
    $bitmapData.Scan0,
    $pixelBytes,
    0,
    $pixelBytes.Length
  )
  $bitmap.UnlockBits($bitmapData)

  $mask = [byte[]]::new($maskSize * $maskSize)
  for ($y = 0; $y -lt $maskSize; $y++) {
    for ($x = 0; $x -lt $maskSize; $x++) {
      $alphaOffset = $y * [Math]::Abs($bitmapData.Stride) + $x * 4 + 3
      if ($pixelBytes[$alphaOffset] -ge 96) {
        $mask[$y * $maskSize + $x] = 1
      }
    }
  }

  $fileName = '{0:X4}.mask' -f [int][char]$character
  [System.IO.File]::WriteAllBytes(
    [System.IO.Path]::Combine($OutputDirectory, $fileName),
    $mask
  )

  $path.Dispose()
  $graphics.Dispose()
  $bitmap.Dispose()
}

$fontFamily.Dispose()
$format.Dispose()
Write-Host "Rendered $($characters.Count) masks to $OutputDirectory"
