# Library Imports - Package names
packages <- c("MASS", "xts", "quantmod", "PortfolioAnalytics",
              "ROI", "ROI.plugin.glpk", "ROI.plugin.quadprog", "tidyverse", 
              "highcharter", "PerformanceAnalytics")

# Install packages if not yet installed - Neat trick from Antoine Soetewey's blog
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

#' import data for each symbol and convert to monthly prices
#' Note conversion to monthly automatically removes NA-VALUES
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

# Exploratory Data Analysis - Correlation Matrix -- Not really need for this analysis however
corr.matrix <- fortify.zoo(etf.rets) %>% ungroup() %>%
            select(-Index) %>%  
            ggpairs()
corr.matrix

###### Portfolio Weighting

#  Portfolio Object Specification - Set up portfolio with objective and constraints(Long Only)

#'  GMV = Global Minimum Variance  ,    EW  = Equally Weight
#'   RO  = Robust Optimization ,   BL  = Black Litterman

port.spec <- portfolio.spec(assets = colnames(etf.rets)) %>%
                add.objective( type="risk", name="StdDev") %>%
                add.objective( type="return", name="mean") %>%
                add.constraint(type="full_investment")  %>%
                add.constraint( type="long_only") %>%
                add.constraint(type = "box", min=0.05, max=0.4)

port.spec.EW <- portfolio.spec(assets = colnames(etf.rets)) %>%
                  add.constraint(type="full_investment")  %>%
                  add.constraint( type="long_only")  %>%
                  add.constraint(type = "box", min=0.05, max=0.4)

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

#' Find Optimal Weights
#' Rf= 0 
GMV.weights<- optimize.portfolio(R = etf.rets,portfolio= port.spec,maxSR=TRUE,
                         optimize_method = "ROI",trace = TRUE,momentFUN='samp.moments')
print("Optimal Markowitz weights")
print(sprintf("%s : %s", etf.symbols,
              scales::percent(as.numeric(GMV.weights$weights))))

#' The following function is obtained from Ross Bennett's (2018) Paper
#' Estimates Covariance Matrix
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


#' The optimizer is not able to find a unique 'tangential' solution, so remove 
#' Max Sharpe ratio objective.
#' This 'loosened' optimization implementation is then focused on maximizing only portfolio
#' return subject to constraints on weights.
BL.weights<-optimize.portfolio(R = etf.rets,portfolio= port.spec,
                               optimize_method = "ROI",trace = TRUE,momentFUN='bl.moments')
print("Black-Litterman Weights")
print(sprintf("%s : %s", etf.symbols,
              scales::percent(as.numeric(BL.weights$weights))))

# Create List of portfolio objects and display their weights in a barplot
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
# rets = returns
rets<-c()

#' Calculate Returns on Out-of-Sample Data - Without Wealth Progression
#' Combine returns of each portfolio in a single data frame
for (i in 1:ncol(weights)) {
  rets.values <- Return.portfolio(etf.rets.test, weights  = weights[,i],
                                  wealth.index = TRUE)
  rets<-cbind(rets,rets.values)
}

# Combine Returns into a single data frame
amount<- 10000
rets.scale<- rets*amount
date <- as.Date("2021-12-31")
init.investment <- xts(matrix(amount, nrow=1,ncol=4), order.by=date)
rets.scale<- rbind(init.investment, rets.scale)
names(rets.scale) <- c("Markowitz","Equal-Weight","Covariance-Robust Optimization",
                       "Black Litterman")

# Performance - Cumulative Returns Chart
cum.rets.chart<-highchart(type='stock')%>%
  hc_title(text='Hypotheical Growth of $10,000 ') %>%
  hc_add_series(round(rets.scale[,1],2), 
                name=colnames(rets.scale)[1])%>%
  hc_add_series(round(rets.scale[,2],2), 
                name=colnames(rets.scale)[2])%>%
  hc_add_series(round(rets.scale[,3],2), 
                name=colnames(rets.scale)[3])%>%
  hc_add_series(round(rets.scale[,4],2), 
                name=colnames(rets.scale)[4])%>%
  hc_navigator(enabled=FALSE)%>%
  hc_scrollbar(enabled=FALSE) %>%
  hc_add_theme(hc_theme_flat())%>%
  hc_exporting(enabled=TRUE)%>%
  hc_legend(enabled=TRUE)
cum.rets.chart

key.measures<- function( Ra, port.weight, bm, scale = 12,Sigma = cov(etf.rets), Rf= 0){
  ret<-     Return.annualized(Ra,scale = scale) *100
  act.ret<-  ActiveReturn(Ra ,Rb= bm,scale = scale) *100
  ex.ret<-   sum(Return.excess(Ra))*100   # Rf= 0 
  #' Note ret = ex.ret since Rf=0
  #' That is, total returns should be roughly equal to excess returns
  #' Where total returns = annualized returns if the investment period < or = 1 year
  st.dev<-  StdDev.annualized(Ra,weights=port.weight,scale = scale) *100
  beta<- CAPM.beta(Ra,bm)
  sharpe<-  SharpeRatio.annualized(Ra, Rf= Rf, scale = scale)
  treynor<- TreynorRatio(Ra, Rb = bm,scale=scale)
  i.ratio<- InformationRatio(Ra, Rb = bm,scale=scale)
  m.square<- Modigliani(Ra, Rb = bm,scale=scale)
  t.error<- TrackingError(Ra, Rb = bm,scale = scale)
  dr <-     FRAPO::dr(port.weight, Sigma ) # Function obtained from FRAPO package
  cr <-     FRAPO::cr(port.weight, Sigma)  # Function obtained from FRAPO package
  res<-     round(rbind(ret,act.ret,ex.ret, st.dev,beta,sharpe, treynor,m.square,i.ratio, t.error, dr,cr),3)
  
  rownames(res) <- c("Annualized Return(%)","Active Return(%)","Excess Return(%)",
                     'Annualized Standard Deviation(%)',"Beta",'Sharpe','Treynor',
                     "M^2 of Modigliani", "Information Ratio","Tracking Error",
                     "Diversification Ratio","Concentration Ratio")
  return(res)
}

#'   Note the portfolio by construction is invested in the overall market.
#'   The Benchmark used here is a Price-Weighted index.
#'   Price- Weighted indices by construction average prices equally.
#'   comp.ind = composite index returns

#'   The motivation for this multi-asset composite index is the
#'   FTSE Multi-Asset Composite Index Series which 
#'   according to FTSE Russell
#'   is designed to measure cross-asset market returns for a range of risk exposures.

etf.prices$comp<-rowMeans(etf.prices)
#' Calculate Returns and remove na value in first row
#' Select Returns from 2022 
comp.ind <- Return.calculate(etf.prices$comp,method='log') %>%
  na.omit()
comp.ind<- comp.ind['2022']

# rets = returns
rets<-c()

#' Calculate Returns on Out-of-Sample Data - Without Wealth Progression
#' Combine returns of each portfolio in a single data frame
for (i in 1:ncol(weights)) {
  rets.values <- Return.portfolio(etf.rets.test, weights  = weights[,i])
  rets<-cbind(rets,rets.values)
}

# ports = portfolios
ports<-c()
# Combine performance measures of each portfolio in a single data frame
for (i in 1:ncol(weights)) {
  measures <- key.measures(port.weight = weights[,i], 
                           Ra = rets[,i], bm =  comp.ind)
  ports<-cbind(ports,measures)
}

#'  Note for the benchmark, all we need are the annualized
#'  Returns, StDev, and Active Returns
#'  Since the benchmark is price-weighted, the weights used in the function is the Equal Weight
bm <- key.measures(port.weight = EW.weights$weights, 
                   Ra = comp.ind, bm =  comp.ind)

keep.rows<- c("Annualized Return(%)","Excess Return(%)",
              'Annualized Standard Deviation(%)')
# Assign NA to those values not in keep.rows
bm[!(bm %in% bm[keep.rows,]),]<- NA 

colnames(all.ports) <- append(names(rets.scale),"Benchmark")

##### Component Contribution to Portfolio Volatility/ Standard Deviation

#   What percentage of portfolio St.Dev did each position account for?

weights <- cbind(GMV.weights$weights,EW.weights$weights,RO.weights$weights,BL.weights$weights)
all.comps.sd<-c()
# Combine component contributions of each portfolio in a single data frame
for (i in 1:ncol(weights)) {
  # Scale by 100 to get percent from decimals like say 0.2 to 20%
  Std <- (StdDev(etf.rets.test, weights = weights[,i], portfolio_method=  
                   "component")$pct_contrib_StdDev)*100
  all.comps.sd<-cbind(all.comps.sd,Std)
}

colnames(all.comps.sd) <- names(rets.scale)
print("Portfolio St.dev Contributions (%)")
all.comps.sd