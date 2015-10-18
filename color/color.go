package color

import (
	"fmt"
	"math"
)

const (
	epsilon float64 = 216.0 / 24389.0
	kappa   float64 = 24389.0 / 27.0
)

// 3 компонента для определения цветов
type Point [3]float64
type Matrix [9]float64

func (p Point) mulMatrix(m *Matrix) (result Point) {
	result[0] = p[0]*m[0] + p[1]*m[1] + p[2]*m[2]
	result[1] = p[0]*m[3] + p[1]*m[4] + p[2]*m[5]
	result[2] = p[0]*m[6] + p[1]*m[7] + p[2]*m[8]
	return
}

func sqr(v float64) float64 {
	return v * v
}

func (rgb RGBColor) String() string {
	return fmt.Sprintf("RGB:{%d %d %d}", rgb.R(), rgb.G(), rgb.B())
}

func (lab LabColor) String() string {
	return fmt.Sprintf("Lab:{%4.7f %4.7f %4.7f}", lab[0], lab[1], lab[2])
}

func (xyz XYZColor) String() string {
	return fmt.Sprintf("XYZ:{%4.7f %4.7f %4.7f}", xyz[0], xyz[1], xyz[2])
}

//Observer = 2°, Illuminant = D65

// Излучатель и матрицы конвертации задают параметры конвертации между цветовыми пространствами RGB и XYZ. По умолчанию установлены на преобразование в пространство sRGB с излучателем D65 и наблюдателем под углов в 2 градуса
var (
	// цветовые координаты излучателя белого света D65 в пространстве XYZ
	WhitePoint = [3]float64{0.95047, 1.0, 1.08883}

	// матрица конвертации из пространства sRGB в пространство XYZ
	ConvertMatrix = Matrix{
		0.412383, 0.357585, 0.18048,
		0.212635, 0.71517, 0.072192,
		0.01933, 0.119195, 0.950528,
	}

	// матрица конвертации из пространства XYZ в пространство sRGB
	InverseMatrix = Matrix{
		3.24103, -1.53741, -0.49862,
		-0.969242, 1.87596, 0.041555,
		0.055632, -0.203979, 1.05698,
	}
)

// коэффициэнт масштабирования RGB
var RGBScale = 255.0

// RGBColor определяет цвет в пространстве RGB, хранящийся в виде 3-х значений в диапазоне от 0 до 1
type RGBColor Point

// R возвращает Red-компоненту RGB-цвета в диапазоне от 0 до RGBScale
func (rgb RGBColor) R() int {
	return int(rgb[0] * RGBScale)
}

// G возвращает Green-компоненту RGB-цвета в диапазоне от 0 до RGBScale
func (rgb RGBColor) G() int {
	return int(rgb[1] * RGBScale)
}

// B возвращает Blie-компоненту RGB-цвета в диапазоне от 0 до RGBScale
func (rgb RGBColor) B() int {
	return int(rgb[2] * RGBScale)
}

// NewRGB создает RGB-цвет на основе 3-х компонент со значениями от 0 до RGBScale
func NewRGB(r, g, b int) (rgb RGBColor) {
	rgb[0] = float64(r) / RGBScale
	rgb[1] = float64(g) / RGBScale
	rgb[2] = float64(b) / RGBScale
	return
}

// RGB конвертирует цвет в пространство RGB
func (rgb RGBColor) RGB() RGBColor {
	return rgb
}

// XYZ конвертирует цвет в пространство XYZ
func (rgb RGBColor) XYZ() XYZColor {
	// linearization
	f := func(n float64) float64 {
		if n > 0.04045 {
			return math.Exp(2.4 * math.Log((n+0.055)/1.055))
		} else {
			return n / 12.92
		}
	}

	p := Point{f(rgb[0]), f(rgb[1]), f(rgb[2])}.mulMatrix(&ConvertMatrix)

	return XYZColor(p)
}

// Lab конвертирует цвет в пространство Lab
func (rgb RGBColor) Lab() LabColor {
	return rgb.XYZ().Lab()
}

type XYZColor Point

// RGB конвертирует цвет в пространство RGB
func (xyz XYZColor) RGB() RGBColor {
	// unlinearization
	f := func(n float64) float64 {
		if n > 0.0031308 {
			return 1.055*math.Exp((1/2.4)*math.Log(n)) - 0.055
		} else {
			return 12.92 * n
		}
	}

	p := Point(xyz).mulMatrix(&InverseMatrix)
	p[0] = f(p[0])
	p[1] = f(p[1])
	p[2] = f(p[2])

	return RGBColor(p)
}

// XYZ конвертирует цвет в пространство XYZ
func (xyz XYZColor) XYZ() XYZColor {
	return xyz
}

// Lab конвертирует цвет в пространство Lab
func (xyz XYZColor) Lab() LabColor {

	f := func(n float64) float64 {
		if n > epsilon {
			return math.Pow(n, 1.0/3.0)
		} else {
			return (kappa*n + 16.0) / 116.0
		}
	}

	x := f(xyz[0] / WhitePoint[0])
	y := f(xyz[1] / WhitePoint[1])
	z := f(xyz[2] / WhitePoint[2])

	return LabColor{
		(116.0 * y) - 16.0,
		500.0 * (x - y),
		200.0 * (y - z),
	}
}

type LabColor Point

func (lab LabColor) L() float64 {
	return lab[0]
}

func (lab LabColor) A() float64 {
	return lab[1]
}

func (lab LabColor) B() float64 {
	return lab[2]
}

// RGB конвертирует цвет в пространство RGB
func (lab LabColor) RGB() RGBColor {
	return lab.XYZ().RGB()
}

// XYZ конвертирует цвет в пространство XYZ
func (lab LabColor) XYZ() XYZColor {

	f := func(n float64) float64 {
		if n3 := n * n * n; n3 > epsilon {
			return n3
		} else {
			return (n*116.0 - 16.0) / kappa
		}
	}

	y := (lab[0] + 16.0) / 116.0
	x := 0.002*lab[1] + y
	z := y - 0.005*lab[2]

	return XYZColor{
		f(x) * WhitePoint[0],
		f(y) * WhitePoint[1],
		f(z) * WhitePoint[2],
	}
}

// Lab конвертирует цвет в пространство Lab
func (lab LabColor) Lab() LabColor {
	return lab
}

func (lab LabColor) LCh() LChColor {
	c := math.Sqrt(sqr(lab[0]) + sqr(lab[2]))
	h := 180.0 * math.Atan2(lab[2], lab[1]) / math.Pi
	if h < 0.0 {
		h += 360.0
	}
	return LChColor{lab[0], c, h}
}

type LChColor Point

// RGB конвертирует цвет в пространство RGB
func (lch LChColor) RGB() RGBColor {
	return lch.Lab().XYZ().RGB()
}

// XYZ конвертирует цвет в пространство XYZ
func (lch LChColor) XYZ() XYZColor {
	return lch.Lab().XYZ()
}

// Lab конвертирует цвет в пространство Lab
func (lch LChColor) Lab() LabColor {
	a := lch[1] * math.Cos(lch[2]*math.Pi/180.0)
	b := lch[1] * math.Sin(lch[2]*math.Pi/180.0)
	return LabColor{lch[0], a, b}
}

type Colorer interface {
	RGB() RGBColor
	XYZ() XYZColor
	Lab() LabColor
}

func ToRGB(color Colorer) RGBColor {
	switch c := color.(type) {
	case XYZColor:
		return c.RGB()
	case LabColor:
		return c.XYZ().RGB()
	case LChColor:
		return c.Lab().XYZ().RGB()
	}
	return color.(RGBColor)
}

type Comparer interface {
	Compare(c1, c2 Colorer) float64
}

type CIE76 struct{}

var DeltaCIE76 CIE76

// DeltaCIE76.Compare computes the CIE76 color difference.
// This is just Euclidean distance in Lab space, and therefore quite fast,
// though it exhibits perceptual uniformity issues especially in the blue and desaturated regions.
func (CIE76) Compare(c1, c2 Colorer) float64 {
	lab1 := c1.Lab()
	lab2 := c2.Lab()
	return math.Sqrt(sqr(lab1[0]-lab2[0]) + sqr(lab1[1]-lab2[1]) + sqr(lab1[2]-lab2[2]))
}

type CIE94 struct {
	// struct for weighting factors for CIE94 ΔE calculation.
	KL, KC, Kh, K1, K2 float64
}

// KLCH94GraphicArts are the weighting factors for CIE94 used for most uses except textiles.
var DeltaCIE94GraphicArts = CIE94{1, 1, 1, 0.045, 0.015}

// KLCH94Textiles are the weighting factors for CIE94 used for textiles.
var DeltaCIE94Textiles = CIE94{2, 1, 1, 0.048, 0.014}

// DeltaCIE94.Compare computes the CIE94 color difference of two L*a*b* colors.
// This is a distance calculation with the addition of weighting factors specified by klch.
func (cie94 CIE94) Compare(c1, c2 Colorer) float64 {
	lab1 := c1.Lab()
	lab2 := c2.Lab()

	dL := sqr(lab1.L() - lab2.L())

	c1ab := math.Sqrt(sqr(lab1.A()) + sqr(lab1.B()))
	c2ab := math.Sqrt(sqr(lab2.A()) + sqr(lab2.B()))

	dC := sqr(c1ab - c2ab)

	dH := sqr(lab1.A()-lab2.A()) + sqr(lab1.B()-lab2.B()) - dC

	sC := 1.0 + cie94.K1*c1ab
	sH := 1.0 + cie94.K2*c1ab

	return math.Sqrt(
		dL/sqr(cie94.KL) +
			dC/sqr(cie94.KC*sC) +
			dH/sqr(cie94.Kh*sH))
}

type CIE2000 struct {
	KL, KC, Kh float64
}

// KLCHDefault is the most commonly used set of weighting parameters for CIEDE2000
var DeltaCIE2000 = CIE2000{1, 1, 1}

// CIE2000 computes the CIEDE2000 delta-E for two L*a*b* space color coordinates
// klch is for configuring the weighting factors, but this almost always should be KLCHDefault
// Note that this implementation will exhibit slightly different behavior around the discontinuities
// of the function (these are grey colors) compared to Java and most C runtimes. The golang atan
// function has different accuracy characteristics compared to most Unix platforms and Java Strict math
func (cie2k CIE2000) Compare(col1, col2 Colorer) float64 {
	lab1 := col1.Lab()
	lab2 := col2.Lab()

	lBarPrime := (lab1.L() + lab2.L()) * 0.5
	c1 := math.Sqrt(lab1.A()*lab1.A() + lab1.B()*lab1.B())
	c2 := math.Sqrt(lab2.A()*lab2.A() + lab2.B()*lab2.B())
	cBar := (c1 + c2) * 0.5

	cBar7 := cBar * cBar * cBar
	cBar7 *= cBar7 * cBar
	g := 0.5 * (1.0 - math.Sqrt(cBar7/(cBar7+6103515625.0))) // 25**7

	a1Prime := (1.0 + g) * lab1.A()
	a2Prime := (1.0 + g) * lab2.A()

	c1Prime := math.Sqrt(a1Prime*a1Prime + lab1.B()*lab1.B())
	c2Prime := math.Sqrt(a2Prime*a2Prime + lab2.B()*lab2.B())

	cBarPrime := (c1Prime + c2Prime) * 0.5

	h1Prime := math.Atan2(lab1.B(), a1Prime)
	if h1Prime < 0 {
		h1Prime += 2 * math.Pi
	}
	h2Prime := math.Atan2(lab2.B(), a2Prime)
	if h2Prime < 0 {
		h2Prime += 2 * math.Pi
	}

	hBarPrime := (h1Prime + h2Prime) * 0.5
	dhPrime := h2Prime - h1Prime
	if math.Abs(dhPrime) > math.Pi {
		hBarPrime += math.Pi
		if h2Prime <= h1Prime {
			dhPrime += 2 * math.Pi
		} else {
			dhPrime -= 2 * math.Pi
		}
	}

	t := 1.0 -
		0.17*math.Cos(hBarPrime-math.Pi/6) +
		0.24*math.Cos(2.0*hBarPrime) +
		0.32*math.Cos(3.0*hBarPrime+math.Pi/30) -
		0.20*math.Cos(4.0*hBarPrime-63.0*math.Pi/180)

	dLPrime := lab2.L() - lab1.L()
	dCPrime := c2Prime - c1Prime
	dHPrime := 2.0 * math.Sqrt(c1Prime*c2Prime) * math.Sin(dhPrime/2.0)

	lBarPrimeM50Sqr := lBarPrime - 50.0
	lBarPrimeM50Sqr *= lBarPrimeM50Sqr
	sL := 1.0 + (0.015*lBarPrimeM50Sqr)/math.Sqrt(20.0+lBarPrimeM50Sqr)
	sC := 1.0 + 0.045*cBarPrime
	sH := 1.0 + 0.015*cBarPrime*t

	hBarPrimeM := (180/math.Pi*hBarPrime - 275.0) / 25.0
	dTheta := math.Pi / 6 * math.Exp(-hBarPrimeM*hBarPrimeM)
	cBarPrime7 := cBarPrime * cBarPrime * cBarPrime
	cBarPrime7 *= cBarPrime7 * cBarPrime
	rC := math.Sqrt(cBarPrime7 / (cBarPrime7 + 6103515625.0))
	rT := -2.0 * rC * math.Sin(2.0*dTheta)

	return math.Sqrt(
		sqr(dLPrime/(cie2k.KL*sL)) +
			sqr(dCPrime/(cie2k.KC*sC)) +
			sqr(dHPrime/(cie2k.Kh*sH)) +
			(dCPrime/(cie2k.KC*sC))*(dHPrime/(cie2k.Kh*sH))*rT)
}

func Approximate(color Colorer, comparer Comparer) (float64, int) {
	best := 0
	bestdist := 10000000.0
	for i, applicant := range XtermLabPalette {
		if dist := comparer.Compare(color, applicant); dist < bestdist {
			best, bestdist = i, dist
		}
	}
	return bestdist, best + 17
}

//color = round(36 * (r * 5) + 6 * (g * 5) + (b * 5) + 16)
func HackApproximate(color Colorer) int {
	c := color.RGB()
	return int(36*(c[0]*5)+6*(c[1]*5)+c[2]*5) + 17
}

var (
	XtermRGBPalette [240]RGBColor
	XtermLabPalette [240]LabColor
)

func init() {
	// calculate xterm 240 color RGB palette

	// calculate 6x6x6 color cube
	cubeLevels := [6]int{0x00, 0x5f, 0x87, 0xaf, 0xd7, 0xff}
	for r := 0; r < 6; r++ {
		for g := 0; g < 6; g++ {
			for b := 0; b < 6; b++ {
				XtermRGBPalette[r*36+g*6+b] = NewRGB(cubeLevels[r], cubeLevels[g], cubeLevels[b])
			}
		}
	}
	// calculate grayscale ramp
	for i := 0; i < 24; i++ {
		XtermRGBPalette[i+216] = NewRGB(0x08+i*0xA, 0x08+i*0xA, 0x08+i*0xA)
	}

	// calculate xterm 240 color Lab palette
	for i, color := range XtermRGBPalette {
		XtermLabPalette[i] = color.Lab()
	}
}
