# Library Imports
library(MASS)
library(xts)
library(quantmod)
library(PortfolioAnalytics)
library(ROI)
library(ROI.plugin.glpk)
library(ROI.plugin.quadprog)
library(tidyverse)
library(highcharter)
library(PerformanceAnalytics)
library(GGally)

# import data for each symbol and convert to monthly prices
# Note conversion to monthly automatically removes NA-VALUES
etf.symbols <- c("SPY", "IWM", "FEZ", "ACWI","AGG","IAGG", "IYR","REET","GSG")
adj.close <- 6  # 6th field is adjusted close
etf.prices <- getSymbols(etf.symbols[1], source="yahoo",
                         auto.assign=FALSE, return.class="xts")[,adj.close]
for (i in 2:length(etf.symbols)) {
  etf.tmp <- getSymbols(etf.symbols[i], source="yahoo",
                        auto.assign=FALSE, return.class="xts")[,adj.close]
  etf.prices <- cbind(etf.prices, etf.tmp)
}
colnames(etf.prices) <- etf.symbols

etf.prices<- suppressWarnings(to.monthly(etf.prices,
                        indexAt='lastof',
                        OHLC=FALSE))

# Calculate Returns and remove na value in first row
etf.rets.data <- Return.calculate(na.omit(etf.prices),method='log') %>%
                na.omit()

# Keep Returns from 2016-2021 for in-sample estimation
etf.rets<- etf.rets.data["2016/2021"]

# Keep Returns from 2022 for out-of-sample testing
etf.rets.test<- etf.rets.data["2022"]

# Exploratory Data Analysis - Correlation Matrix
corr.matrix <- fortify.zoo(etf.rets) %>% ungroup() %>%
            select(-Index) %>%  
            ggpairs()
corr.matrix

###### Portfolio Weighting

#  Portfolio Object Specification - Set up portfolio with objective and constraints(Long Only)

#  GMV = Global Minimum Variance  ,    EW  = Equally Weight
#   RO  = Robust Optimization ,   BL  = Black Litterman

port.spec <- portfolio.spec(assets = colnames(etf.rets)) %>%
                add.objective( type="risk", name="StdDev") %>%
                add.objective( type="return", name="mean") %>%
                add.constraint(type="full_investment")  %>%
                add.constraint( type="long_only")  %>%
                add.constraint(type = "box", min=0.05, max=0.5)

port.spec.EW <- portfolio.spec(assets = colnames(etf.rets)) %>%
                  add.constraint(type="full_investment")  %>%
                  add.constraint( type="long_only")  %>%
                  add.constraint(type = "box", min=0.05, max=0.5)

# Set Portfolio Moments for Optimization Functions
samp.moments<-  function(R){
    out<-set.portfolio.moments(R = R, portfolio=port.spec, method= 'sample')
    return(out)
}


# Views are generated from a VAR-2 model from the MTS package
# R is an xts return object
bl.moments<- function(R){
  data<- fortify.zoo(R)[,-1]
  m1<- MTS::VAR(data,p=2,output = F) # From MTS package
  pred<- MTS::VARpred(m1,1,output = F) # 1 period ahead forecast
  pred<- pred$pred   
  n.assets<-ncol(R)
  pick<-diag(n.assets)
  out<- black.litterman(R= R,P = pick,Views = pred)
  out$mu<- out$BLMu
  out$sigma<- out$BLSigma
  return(out)
}

# Find Optimal Weights
# Rf= 0 
GMV.weights<- optimize.portfolio(R = etf.rets,portfolio= port.spec,maxSR=TRUE,
                         optimize_method = "ROI",trace = TRUE,momentFUN='samp.moments')
print("Optimal Markowitz weights")
print(sprintf("%s : %s", etf.symbols,
              scales::percent(as.numeric(GMV.weights$weights))))

# The following function is obtained from Ross Bennett's (2018) Paper
# Estimates Covariance Matrix
sigma.robust<- function(R){
  out<- list()
  set.seed(1234)
  out$sigma<- cov.rob(R,method='mcd')$cov
  return(out)
}

RO.weights <-optimize.portfolio(R = etf.rets,portfolio= port.spec,maxSR=TRUE,
                                 optimize_method = "ROI",trace = TRUE,momentFUN='sigma.robust')

print("Covariance-Robust Optimization Weights")
print(sprintf("%s : %s", etf.symbols,
              scales::percent(as.numeric(RO.weights$weights))))


EW.weights <- equal.weight(etf.rets, portfolio= port.spec.EW)
print("Equal weights")
print(sprintf("%s : %s", etf.symbols,
              scales::percent(as.numeric(EW.weights$weights))))

BL.weights<-optimize.portfolio(R = etf.rets,portfolio= port.spec,
                               optimize_method = "ROI",trace = TRUE,momentFUN='bl.moments')
print("Black-Litterman Weights")
print(sprintf("%s : %s", etf.symbols,
              scales::percent(as.numeric(BL.weights$weights))))

# Create List of portfolios objects and display their weights in a barplot
port.list <- list(GMV.weights,EW.weights, RO.weights, BL.weights)
names(port.list) <- c("Markowitz","Equal-Weight", 'Covariance Robust Optimization', 
                      'Black-Litterman')
chart.Weights(combine.optimizations(port.list), plot.type = "barplot", 
              main="Weights stacked by Asset (Long-only)")

# Put Weights in a data frame(Scale by 100 to get it in percentages)
weights<- cbind(GMV.weights$weights,EW.weights$weights,
                RO.weights$weights, BL.weights$weights)
weights<- round(weights,4)*100
colnames(weights)<- names(port.list)
print("Weights in (%)")
weights

###### Portfolio Performance Evaluation

# Calculate Returns on Out-of-Sample Data - With Wealth Progression
GMV.rets<- Return.portfolio(etf.rets.test, weights  = unlist(GMV.weights['weights']),
                            wealth.index = TRUE)%>%
          'colnames<-' ('Returns')

EW.rets<-  Return.portfolio(etf.rets.test, weights  = unlist(EW.weights['weights']),
                            wealth.index = TRUE)%>%
         'colnames<-' ('Returns')

RO.rets<-  Return.portfolio(etf.rets.test, weights  = unlist(RO.weights['weights']),
                            wealth.index = TRUE)%>%
        'colnames<-' ('Returns')

BL.rets<-  Return.portfolio(etf.rets.test, weights  = unlist(BL.weights['weights']),
                            wealth.index = TRUE)%>%
        'colnames<-' ('Returns')

# Combine Returns into a single data frame
amount<- 10000
combined.rets<- cbind(GMV.rets$Returns,EW.rets$Returns,RO.rets$Returns,BL.rets$Returns) *amount
date <- as.Date("2021-12-31")
init.investment <- xts(matrix(amount, nrow=1,ncol=4), order.by=date)
combined.rets<- rbind(init.investment, combined.rets)
names(combined.rets) <- c("Markowitz","Equal-Weight","Covariance-Robust Optimization",
                          "Black Litterman")

# Performance - Cumulative Returns Chart
cum.rets.chart<-highchart(type='stock')%>%
                hc_title(text='Hypotheical Growth of $10,000 ') %>%
                hc_add_series(combined.rets[,1], 
                              name=colnames(combined.rets)[1])%>%
                hc_add_series(combined.rets[,2], 
                               name=colnames(combined.rets)[2])%>%
                hc_add_series(combined.rets[,3], 
                              name=colnames(combined.rets)[3])%>%
                hc_add_series(combined.rets[,4], 
                              name=colnames(combined.rets)[4])%>%
                hc_navigator(enabled=FALSE)%>%
                hc_scrollbar(enabled=FALSE) %>%
                hc_add_theme(hc_theme_flat())%>%
                hc_exporting(enabled=TRUE)%>%
                hc_legend(enabled=TRUE)
cum.rets.chart


# Key Performance Measures
# dr  = diversification ratio,  cr  = concentration ratio
# Ra = portfolio returns , bm = benchmark
# port.weight = portfolio weights

key.measures<- function( Ra, port.weight, bm, scale = 12,Sigma = cov(etf.rets), Rf= 0){
  ret<-     Return.annualized(Ra,scale = scale) *100
  act.ret<-  ActiveReturn(Ra ,Rb= bm,scale = scale) *100
  ex.ret<-   sum(Return.excess(Ra))*100   # Rf= 0 
  st.dev<-  StdDev.annualized(Ra,weights=port.weight,scale = scale) *100
  sharpe<-  SharpeRatio.annualized(Ra, Rf= Rf, scale = scale)
  treynor<- TreynorRatio(Ra, Rb = bm,scale=scale)
  i.ratio<- InformationRatio(Ra, Rb = bm,scale=scale)
  m.square<- Modigliani(Ra, Rb = bm,scale=scale)
  t.error<- TrackingError(Ra, Rb = bm,scale = scale)
  dr <-     FRAPO::dr(port.weight, Sigma ) # Function obtained from FRAPO package
  cr <-     FRAPO::cr(port.weight, Sigma)  # Function obtained from FRAPO package
  res<-     rbind(ret,act.ret,ex.ret, st.dev,sharpe, treynor,i.ratio,m.square, t.error, dr,cr) 
  rownames(res) <- c("Return(%)","Active Return(%)","Excess Return(%)",
                     'Standard Deviation(%)','Sharpe','Treynor',
                  "M^2 of Modigliani", "Information Ratio","Tracking Error",
                  "Diversification Ratio","Concentration Ratio")
  return(res)
}

#   Note the portfolio by construction is invested in the overall market.
#   The Benchmark used here is a Price-Weighted index.
#   Price- Weighted indices by construction average prices equally.
#   comp.ind = composite index returns

#   The motivation for this multi-asset composite index is the
#   FTSE Multi-Asset Composite Index Series which 
#   according to FTSE Russell
#   is designed to measure cross-asset market returns for a range of risk exposures.

etf.prices$comp<-rowMeans(etf.prices)
# Calculate Returns and remove na value in first row
# Select Returns from 2022 
comp.ind <- Return.calculate(etf.prices$comp,method='log') %>%
            na.omit()
comp.ind<- comp.ind['2022']

# Calculate Returns on Out-of-Sample Data - Without Wealth Progression
GMV.rets<- Return.portfolio(etf.rets.test, weights  = GMV.weights$weights)

EW.rets<-  Return.portfolio(etf.rets.test, weights  = EW.weights$weights)

RO.rets<-  Return.portfolio(etf.rets.test, weights  = RO.weights$weights)

BL.rets<-  Return.portfolio(etf.rets.test, weights  = BL.weights$weights)


gmv<-key.measures(port.weight = GMV.weights$weights, 
                  Ra = GMV.rets[,1], bm =  comp.ind)
ew<-key.measures(port.weight = EW.weights$weights, 
                 Ra = EW.rets[,1], bm =  comp.ind)
ro<-key.measures(port.weight = RO.weights$weights, 
                 Ra = RO.rets[,1], bm =  comp.ind)
bl<-key.measures(port.weight = BL.weights$weights, 
                 Ra = BL.rets[,1], bm =  comp.ind)

#  Note for the benchmark, all we need are the annualized
#  Returns, StDev, and Sharpe Ratios
#  Since the benchmark is price-weighted, the weights used in the function is the Equal Weight
bm <- key.measures(port.weight = EW.weights$weights, 
                  Ra = comp.ind, bm =  comp.ind)

bm[c(2,3,6,7,8,9,10,11)]<- NA
all.ports <- cbind(gmv,ew,ro,bl,bm)

colnames(all.ports) <- c("Markowitz","Equal-Weight","Covariance-Robust Optimization",
                         "Black Litterman","Benchmark")
print("Key Portfolio Performance Measures - Annualized")
all.ports

##### Component Contribution to Portfolio Volatility/ Standard Deviation

#   What percentage of portfolio St.Dev did each position account for?

gmv<- StdDev(etf.rets.test, weights =GMV.weights$weights, portfolio_method = "component")
gmv<- round(gmv$pct_contrib_StdDev,3)

ew<- StdDev(etf.rets.test, weights = EW.weights$weights, portfolio_method = "component")
ew<- round(ew$pct_contrib_StdDev,3)

ro<- StdDev(etf.rets.test, weights =RO.weights$weights, portfolio_method = "component")
ro<- round(ro$pct_contrib_StdDev,3)

bl<- StdDev(etf.rets.test, weights =BL.weights$weights, portfolio_method = "component")
bl<- round(bl$pct_contrib_StdDev,3)

# Combine component contributions of each portfolio in a single data frame
all.comps.sd <- cbind(gmv,ew,ro,bl)*100

colnames(all.comps.sd) <- c("Markowitz","Equal-Weight","Covariance-Robust Optimization",
                         "Black Litterman")
print("Portfolio St.dev Contributions (%)")
all.comps.sd