# Backtesting-Expected-Shortfall

## Abstract

With the Fundamental Review of the Trading Book, the Basel Committee on Banking Supervision proposed to replace the Value-at-Risk as the risk measure with the Expected Shortfall in order to calculate the regulatory capital. Since then, both risk measures have to be calculated but model validation and backtesting continues to be based on the Value-at-Risk. In the last years, an academic and practical debate ensued on whether and how the Expected Shortfall is backtestable. One central question is whether the property of elicitability, which is not fulfilled by the Expected Shortfall, is necessary to construct suitable backtests. Over the years, there have been a number of different backtest proposals. In this thesis, six selected backtests are examined in an exemplary manner and integrated into a theoretical framework. One backtest, namely a proposal of Acerbi and Szekely (2017), stands out. The reasons for this are, on the one hand, its ease of implementing. On the other hand, the small sensitivity to Value-at-Risk predictions, if only the Expected Shortfall should be backtested. Based on this test, a traffic light approach similar to the previous Value-at-Risk backtesting framework can be recommended.

## Keywords

Value-at-Risk, Expected Shortfall, Elicitability, Backtesting,
Model Validation, Regulatory Capital

## Notes

First, open "Libraries and functions", install/load all packages and read in the functions there.
Second, perform analysis.R

Note that not all written functions are in the first script, since some examples are based on VaR & ES Setups which one has to select first. One will see this below.
The structure of the analysis code is the same as in the thesis. Namely:

- Significance (IID)
  - A first impression      
  - A fixed number of VaR exceedances 
- Power (IID)
   - A first impression 
   - Underestimated SD 
   - A fixed number of VaR exceedances 
- Which role does the VaR play? (IID)
   - Correct ES but wrong VaR prediction 
   - Underestimated ES but correct VaR prediction 
   - An example of deliberate deception 
- GARCH, the Non-IID-Setup
   - Another example of underestimated SD 
   - The analysis similar to Du and Escanciano
