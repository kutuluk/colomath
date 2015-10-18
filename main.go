package main

import (
	"fmt"

	"github.com/kutuluk/color/color"
)

func main() {

	//	col1 := color.NewRGB(50, 120, 0)
	//	col2 := color.NewRGB(50, 0, 0)

	col1 := color.LabColor{50, 120, 0}
	//col2 := color.LabColor{50, 0, 0}
	col2 := color.ToRGB(col1)
	col3 := color.LabColor{50, 120, 0}.RGB()

	distance76 := color.DeltaCIE76.Compare(col1, col2)
	distance94 := color.DeltaCIE94GraphicArts.Compare(col1, col2)
	distance2000 := color.DeltaCIE2000.Compare(col1, col2)

	fmt.Printf("col1 %v\n", col1)
	fmt.Printf("col2 %v\n", col2)
	fmt.Printf("col3 %v\n", col3)
	fmt.Println("")
	fmt.Println(col1.RGB())
	fmt.Println(col1.XYZ())
	fmt.Println(col1.Lab())
	fmt.Println("")
	fmt.Printf("d76: %f (=120)\n", distance76)
	fmt.Printf("d94: %f (=18.75)\n", distance94)
	fmt.Printf("d2000: %f (=32.44)\n", distance2000)
}
