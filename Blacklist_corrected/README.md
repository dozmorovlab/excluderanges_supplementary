# Blacklist Annotation Correction

This update refactors the region annotation logic by replacing the integer-based hitCode flag with two separate boolean flags, improving clarity and precision.

## Changes Overview

1.	Initialization of Boolean Flags
    - The single integer `hitCode` is replaced with two boolean variables:
        - `highSig` (for high-signal regions)
        - `lowMapp` (for low-mappability regions)
    - Updated Initialization:
        ```diff
        429c429,430
        - 	int hitCode = 0;
        ---
        + 	bool highSig = false;
        + 	bool lowMapp = false;
        ```

2.	Simplified Assignments
    - The use of bitwise OR (`|`) for setting flags is replaced with straightforward boolean assignments based on conditions.
    - Updated Assignments:
        ```diff
        459c460
        - 				hitCode = hitCode | 1;
        ---
        + 				highSig = true;
        465c466
        - 				hitCode = hitCode | 2;
        ---
        + 				lowMapp = true;
        ```

3.	Enhanced Conditional Logic for Region Classification
    - Instead of evaluating the hitCode integer, the new approach uses the boolean flags to classify regions.
        -	If both flags are `true`, the new annotation is "High Signal, Low Mappability".
        -	If neither flag is `true`, the annotation is left blank.
    - Updates to Conditional Logic:
        ```diff
        476c477,479
        - 					if(hitCode == 2) {
        ---
        + 					if (highSig && lowMapp) {
        + 						regionClass = "High Signal Region, Low Mappability";
        + 					} else if (lowMapp) {
        478c481
        - 					} else {
        ---
        + 					} else if (highSig) {
        479a483,484
        + 					} else {
        + 						regionClass = "";
        483c488,489
        - 					hitCode = 0;
        ---
        + 					highSig = false;
        + 					lowMapp = false;
        ```

4.	Last Region Handling
    - The same logic is applied for the conditional block that ouputs the last output region.
    - Revised Block:

        ```diff
        493c499,501
        - 		if(hitCode == 2) {
        ---
        + 		if (highSig && lowMapp) {
        + 			regionClass = "High Signal Region, Low Mappability";
        + 		} else if (lowMapp) {
        495c503
        - 		} else {
        ---
        + 		} else if (highSig) {
        496a505,506
        + 		} else {
        + 			regionClass = "";
        ```