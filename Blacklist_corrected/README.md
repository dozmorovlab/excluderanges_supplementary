```diff
429c429,430
< 	int hitCode = 0;
---
> 	bool highSig = false;
> 	bool lowMapp = false;
459c460
< 				hitCode = hitCode | 1;
---
> 				highSig = true;
465c466
< 				hitCode = hitCode | 2;
---
> 				lowMapp = true;
476c477,479
< 					if(hitCode == 2) {
---
> 					if (highSig && lowMapp) {
> 						regionClass = "High Signal Region, Low Mappability";
> 					} else if (lowMapp) {
478c481
< 					} else {
---
> 					} else if (highSig) {
479a483,484
> 					} else {
> 						regionClass = "";
483c488,489
< 					hitCode = 0;
---
> 					highSig = false;
> 					lowMapp = false;
493c499,501
< 		if(hitCode == 2) {
---
> 		if (highSig && lowMapp) {
> 			regionClass = "High Signal Region, Low Mappability";
> 		} else if (lowMapp) {
495c503
< 		} else {
---
> 		} else if (highSig) {
496a505,506
> 		} else {
> 			regionClass = "";
```