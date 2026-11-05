## makeFig2.R -- Fig2 v9 (Fig2-only redo).  Panel j replaced: the two top sub-panels
## are now Monocle3-style pseudotime trajectory plots on the host UMAP (replacing the
## old DC-plane trajectory2 images); bottom row keeps the 7-bin state-along panels.
## Global compactness: ink-bbox crop-fill, plot.margin = unit(1,'mm'), target max
## blank band <=4mm; canvas 183mm wide, content-adaptive height <=230mm.
options(warn = 1)
source("analysis/figs_final/_v9/R_scripts/common.R")
suppressMessages({library(dplyr); library(RColorBrewer); library(scales)})
theme_set(base_theme())

FIGS <- "analysis/L1_scRNA/results/figs"
PT   <- "analysis/L1_scRNA/results/pseudotime"
V9   <- "analysis/figs_final/_v9"

## =================================================================================
##  PART A -- regenerate the two Monocle3-style trajectory source PNGs (120x90mm,
##  600dpi, white).  Cells coloured by cluster (Set2); black main path = loess curve
##  through dpt-ordered cluster medians; root = min-decouple cluster (black triangle);
##  branch clusters (off the main curve) joined by black dashed connector + cluster-id
##  label; node dots at each median.  Footnote as caption.
## =================================================================================
dist2poly <- function(x0,y0,cx,cy){
  n <- length(cx); if (n < 2) return(Inf)
  seg <- function(a,b){
    dx<-b[1]-a[1]; dy<-b[2]-a[2]; L2<-dx*dx+dy*dy
    if(L2==0) return(sqrt((x0-a[1])^2+(y0-a[2])^2))
    t<-max(0,min(1,((x0-a[1])*dx+(y0-a[2])*dy)/L2))
    sqrt((x0-(a[1]+t*dx))^2+(y0-(a[2]+t*dy))^2)
  }
  pts<-cbind(cx,cy)
  vapply(seq_len(nrow(pts)-1),function(i) seg(pts[i,],pts[i+1,]),numeric(1))
}
nearest_on_curve <- function(x0,y0,cx,cy){
  n<-length(cx); best<-Inf; bi<-1; bt<-0
  for(i in seq_len(n-1)){
    a<-c(cx[i],cy[i]); b<-c(cx[i+1],cy[i+1]); dx<-b[1]-a[1]; dy<-b[2]-a[2]; L2<-dx*dx+dy*dy
    if(L2==0) next
    t<-max(0,min(1,((x0-a[1])*dx+(y0-a[2])*dy)/L2))
    px<-a[1]+t*dx; py<-a[2]+t*dy; dd<-sqrt((x0-px)^2+(y0-py)^2)
    if(dd<best){best<-dd;bi<-i;bt<-t}
  }
  c(px=cx[bi]+bt*(cx[bi+1]-cx[bi]), py=cy[bi]+bt*(cy[bi+1]-cy[bi]), d=best)
}
fit_curve <- function(dat, span=0.75){
  sp <- seq(min(dat$dpt), max(dat$dpt), length.out=220)
  lx <- loess(u1~dpt, data=dat, span=span); ly <- loess(u2~dpt, data=dat, span=span)
  data.frame(x=predict(lx,newdata=data.frame(dpt=sp)), y=predict(ly,newdata=data.frame(dpt=sp)))
}

traj_plot <- function(host, root_cluster, title, min_thresh, thresh_q=0.75){
  d <- read.csv(file.path(PT, sprintf("PT_%s_min_decouple_cells.csv", host)))
  c <- d |> group_by(cluster) |>
    summarise(u1=median(umap1), u2=median(umap2), dpt=median(dpt), n=n(), .groups="drop") |>
    arrange(dpt)
  c$root <- c$cluster == root_cluster
  ncl <- nrow(c)
  cv  <- fit_curve(c)
  d0  <- vapply(seq_len(nrow(c)), function(i) min(dist2poly(c$u1[i],c$u2[i], cv$x, cv$y)), numeric(1))
  thr <- pmax(min_thresh, quantile(d0, thresh_q))
  br  <- (d0 > thr) & (!c$root)
  cv2 <- fit_curve(c[!br,])                      # main path curve (backbone)
  if (any(br)) {
    conn <- lapply(which(br), function(i) nearest_on_curve(c$u1[i], c$u2[i], cv2$x, cv2$y))
    conn <- do.call(rbind, conn)
    brdf <- data.frame(x0=c$u1[br], y0=c$u2[br], px=conn[,1], py=conn[,2], grp=c$cluster[br])
  } else {
    brdf <- data.frame(x0=numeric(0), y0=numeric(0), px=numeric(0), py=numeric(0), grp=character(0))
  }
  set2 <- RColorBrewer::brewer.pal(8,"Set2")
  pal  <- colorRampPalette(set2)(ncl)
  dd   <- data.frame(umap1=d$umap1, umap2=d$umap2, cluster=factor(d$cluster))
  alpha <- if (nrow(dd) <= 20000) 0.6 else 0.35
  xr <- range(c$u1); yr <- range(c$u2)
  xlim <- c(xr[1]-0.10*diff(xr), xr[2]+0.10*diff(xr))
  ylim <- c(yr[1]-0.12*diff(yr), yr[2]+0.12*diff(yr))
  barcode <- if (host=="fib") "CCCAGTTAGCCAGTAG-50" else "AAATGCCAGGGTATCG-6"
  p <- ggplot() +
    geom_point(data=dd, aes(umap1,umap2,colour=cluster), size=0.7, alpha=alpha) +
    scale_colour_manual(values=pal) +
    geom_path(data=cv2, aes(x=x,y=y), colour="black", linewidth=1.5, lineend="round") +
    geom_segment(data=brdf, aes(x=x0,y=y0,xend=px,yend=py), colour="black",
                 linewidth=0.45, linetype="dashed", alpha=0.85) +
    geom_point(data=c, aes(u1,u2), size=1.8, colour="black") +
    geom_point(data=c[c$root,], aes(u1,u2), shape=17, size=3, colour="black") +
    geom_text(data=cbind(c, br=br)[br,], aes(u1,u2,label=cluster), size=2.4,
              colour="#880000", hjust=-0.35, vjust=-0.35, family="Arial") +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    labs(title=title, x="UMAP1", y="UMAP2",
         caption=paste0("pseudotime on host UMAP; root = min-decouple cluster (cluster ",
               root_cluster, "; barcode ", barcode, "); DC1 limited \u2192 embedded UMAP; ",
               "branch structure from cluster medians")) +
    theme_classic(base_size=6.5, base_family="Arial") +
    theme(plot.title=element_text(size=7,face="bold",hjust=0),
          plot.caption=element_text(size=6, face="plain", colour="#444444", hjust=0),
          axis.title=element_text(size=6.2), axis.text=element_text(size=6.0),
          legend.position="none", plot.margin=margin(2,2,2,2))
  list(p=p, br=c$cluster[br], thr=thr, ncl=ncl, dist=d0)
}

fx <- traj_plot("fib", 14, "Fibroblasts (n=17,908; state 2.38%)", min_thresh=0.25)
px <- traj_plot("plasma", 4, "Plasma cells (n=2,397; state 2.92%)", min_thresh=0.10)
ggsave(file.path(V9,"PT_fib_umap_traj.png"),   fx$p, width=120/25.4, height=90/25.4, dpi=600, bg="white")
ggsave(file.path(V9,"PT_plasma_umap_traj.png"), px$p, width=120/25.4, height=90/25.4, dpi=600, bg="white")
cat(sprintf("traj: fib branches=%s (thr=%.3f); plasma branches=%s (thr=%.3f)\n",
            paste(fx$br,collapse=","), fx$thr, paste(px$br,collapse=","), px$thr))

## ---------------- compact state-along panels (120x70mm, 600dpi, no title block) ----
wilson <- function(k, n, z=1.96){
  if (n==0) return(c(lo=NA, hi=NA))
  p <- k/n; d <- 1 + z^2/n
  ctr <- (p + z^2/(2*n))/d
  hw  <- z*sqrt((p*(1-p) + z^2/(4*n))/n)/d
  c(lo=ctr-hw, hi=ctr+hw)
}
state_along_c_plot <- function(host){
  d <- read.csv(file.path(PT, sprintf("PT_%s_min_decouple_cells.csv", host)))
  d$state <- as.logical(as.character(d$state))
  brk <- quantile(d$dpt, seq(0,1,length.out=8))
  d$bin <- cut(d$dpt, breaks=brk, include.lowest=TRUE, labels=FALSE)
  g <- d |> group_by(bin) |> summarise(n=n(), state=sum(state), dpt_mid=median(dpt), .groups="drop")
  glob <- mean(d$state)
  w <- t(vapply(seq_len(nrow(g)), function(i) wilson(g$state[i], g$n[i]), numeric(2)))
  g$fr <- g$state/g$n; g$lo <- w[,1]; g$hi <- w[,2]; g$sig <- g$lo > glob
  ymin <- max(0, min(g$lo)*100*0.9); ymax <- max(g$hi)*100*1.05
  xr <- range(g$dpt_mid); padx <- 0.06*diff(xr)
  xlim <- c(xr[1]-padx, xr[2]+padx)
  ggplot(g, aes(x=dpt_mid, y=fr*100)) +
    geom_hline(yintercept=glob*100, colour="#555555", linewidth=0.4, linetype=2) +
    geom_errorbar(aes(ymin=lo*100, ymax=hi*100), width=0.012, colour="#3C5488", linewidth=0.5) +
    geom_line(colour="#3C5488", linewidth=0.8) +
    geom_point(size=1.8, colour="#3C5488") +
    geom_point(data=subset(g, sig), aes(x=dpt_mid, y=hi*100), shape=8, size=2.6, colour="#E64B35") +
    scale_x_continuous(limits=xlim, breaks=pretty(xlim, 4),
                       labels=function(x) sprintf("%.2f", x)) +
    scale_y_continuous(limits=c(ymin, ymax), breaks=pretty(c(ymin, ymax)), labels=function(x) sprintf("%.1f", x)) +
    labs(x="diffusion pseudotime (dpt)", y="state fraction (%)") +
    annotate("text", x=xlim[2], y=ymax, hjust=1, vjust=1.2, size=2.15, colour="#333333",
             family="Arial", label="red *= sig bin; bars = Wilson 95% CI") +
    theme_classic(base_size=6.5, base_family="Arial") +
    theme(axis.title=element_text(size=6.2), axis.text=element_text(size=6.0),
          axis.line=element_line(linewidth=0.35),
          plot.margin=margin(1.5,1.5,1.0,1.5,"pt"),
          plot.background=element_rect(fill="white", colour=NA),
          panel.background=element_rect(fill="white", colour=NA))
}
sc_f <- state_along_c_plot("fib");   ggsave(file.path(V9,"PT_fib_state_along_c.png"),
                                            sc_f, width=120/25.4, height=70/25.4, dpi=600, bg="white")
sc_p <- state_along_c_plot("plasma");ggsave(file.path(V9,"PT_plasma_state_along_c.png"),
                                            sc_p, width=120/25.4, height=70/25.4, dpi=600, bg="white")
cat("compact state_along_c written\n")


## =================================================================================
##  PART B -- panels a-i.  a = UMAP hero (raster, kept); b-g = R vector redraws;
##  h,i = R vector (existing).  Text >=6.5pt.
## =================================================================================
SC <- read.csv("analysis/L1_scRNA/results/module_scores/scores_GSE266577.csv")
SC$state <- as.logical(as.character(SC$state_high))
RAR <- read.csv(file.path(FIGS, "F2d_rarity.csv"))
WW  <- SZ <- NULL

## ---- panel b: quadrant (M2+M3 z mean vs M1z; red = state) ----
dq <- data.frame(m23=(SC$M2z+SC$M3z)/2, m1=SC$M1z, state=SC$state)
set.seed(1); dq <- dq[sample(nrow(dq), min(35000, nrow(dq))), ]
p_b <- ggplot(dq[!dq$state,], aes(x=m23, y=m1)) +
  geom_point(size=0.25, alpha=0.20, colour="grey55") +
  geom_point(data=dq[dq$state,], size=0.6, alpha=0.75, colour=SIG) +
  geom_hline(yintercept=0, linewidth=0.3, colour="grey70", linetype=2) +
  geom_vline(xintercept=0, linewidth=0.3, colour="grey70", linetype=2) +
  labs(title="Quadrant enrichment (M1z vs M2/M3z)", x="M2+M3 z (mean)", y="M1z",
       caption="red = cell state (n=1,751)") +
  theme_classic(base_size=6.5, base_family="Arial") +
  theme(plot.title=element_text(size=8,face="bold",hjust=0), plot.caption=element_text(size=6.5,colour="#333333",hjust=0),
        axis.title=element_text(size=6.5), axis.text=element_text(size=6.5), plot.margin=margin(2,2,2,2)) +
  labs(tag="b")

## ---- panel c: rarity (5 source state% + binom CI) ----
z <- 1.96
RAR$lo <- (RAR$state_n + z^2/2)/(RAR$n+z^2) - z*sqrt((RAR$state_n*(1-RAR$state_n/RAR$n)+z^2/4)/RAR$n)/(RAR$n+z^2)
RAR$hi <- (RAR$state_n + z^2/2)/(RAR$n+z^2) + z*sqrt((RAR$state_n*(1-RAR$state_n/RAR$n)+z^2/4)/RAR$n)/(RAR$n+z^2)
RAR$lo <- pmax(0,RAR$lo)
p_c <- ggplot(RAR, aes(x=pct, y=factor(source, levels=rev(source)))) +
  geom_errorbarh(aes(xmin=lo*100, xmax=hi*100), height=0.2, colour=BLU, linewidth=0.5) +
  geom_point(size=2.4, colour=BLU) +
  scale_x_log10(breaks=c(1,1.5,2,2.5), labels=function(x) sprintf("%.2f",x)) +
  labs(title="Rare cell state prevalence by source", x="state prevalence (%) [log]", y=NULL) +
  theme_classic(base_size=6.5, base_family="Arial") +
  theme(plot.title=element_text(size=8,face="bold",hjust=0), axis.title.x=element_text(size=6.5),
        axis.text.y=element_text(size=6.5), axis.text.x=element_text(size=6.5), plot.margin=margin(2,2,2,2), legend.position="none") +
  labs(tag="c")

## ---- panel d: OR forest (state by lineage, 95% CI) ----
lines <- c("Epithelial cells","Mast cells","Endothelial cells","NK cells","Fibroblasts","Plasma cells")
drow <- lapply(lines, function(lg){
  tab <- table(state=SC$state, lineage=(SC$cell_type==lg))
  a2<-tab[2,2];b2<-tab[2,1];c3<-tab[1,2];d3<-tab[1,1]
  or<-(a2*d3)/(b2*c3); se<-sqrt(1/a2+1/b2+1/c3+1/d3)
  c(or=or, lo=exp(log(or)-1.96*se), hi=exp(log(or)+1.96*se))
})
dor <- as.data.frame(do.call(rbind,drow)); dor$lineage<-lines
p_d <- ggplot(dor, aes(x=or, y=factor(lineage, levels=lines))) +
  geom_vline(xintercept=1, linewidth=0.35, colour="grey55", linetype=2) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.2, colour=BLU, linewidth=0.5) +
  geom_point(size=2.4, colour=BLU) +
  scale_x_log10(breaks=c(0.15,0.3,1,3), labels=function(x) sprintf("%.2f",x)) +
  labs(title="State enrichment by lineage (OR)", x="odds ratio (log)", y=NULL) +
  theme_classic(base_size=6.5, base_family="Arial") +
  theme(plot.title=element_text(size=8,face="bold",hjust=0), axis.title.x=element_text(size=6.5),
        axis.text.y=element_text(size=6.5), axis.text.x=element_text(size=6.5), plot.margin=margin(2,2,2,2), legend.position="none") +
  labs(tag="d")

## ---- panel e: 3 violins M1/M2/M3 (state vs rest) + MWU p ----
md <- data.frame(mod=rep(factor(c("M1","M2","M3"),levels=c("M1","M2","M3")), each=nrow(SC)),
                 state=rep(SC$state,3), val=c(SC$M1z,SC$M2z,SC$M3z))
pvs <- vapply(c("M1z","M2z","M3z"), function(m) wilcox.test(SC[[m]][SC$state],SC[[m]][!SC$state])$p.value, numeric(1))
p_e <- ggplot(md, aes(x=mod, y=val, fill=state, group=interaction(mod,state))) +
  geom_violin(position=position_dodge(0.7), width=0.65, scale="area", alpha=0.8, colour="grey30", linewidth=0.2) +
  scale_fill_manual(values=c("TRUE"=SIG,"FALSE"=GRA), name=NULL, labels=c("rest","state")) +
  labs(title="Module z by state (M1/M2/M3)", x=NULL, y="z-score",
       caption="MWU p<1e-300 each (state n=1,751)") +
  theme_classic(base_size=6.5, base_family="Arial") +
  theme(plot.title=element_text(size=8,face="bold",hjust=0), plot.caption=element_text(size=6.5,colour="#333333",hjust=0),
        axis.title.y=element_text(size=6.5), axis.text=element_text(size=6.5),
        legend.position="bottom", legend.text=element_text(size=6.5), legend.key.size=unit(0.22,"cm"),
        plot.margin=margin(2,2,2,2)) +
  labs(tag="e")

## ---- panel f: composition stacked (all vs state), 8 major categories ----
ct <- as.character(SC$cell_type)
mkcat <- function(cc){
  o <- rep("Other", length(cc))
  o[grepl("T cells|Tem|Tcm|Regulatory|Type 17", cc)] <- "T cells"
  o[grepl("Macrophage|monocyte|DC|pDC", cc)]        <- "Myeloid"
  o[grepl("NK", cc)]                                <- "NK cells"
  o[grepl("B cells", cc)]                           <- "B cells"
  o[grepl("Plasma", cc)]                            <- "Plasma cells"
  o[grepl("Epithelial", cc)]                        <- "Epithelial cells"
  o[grepl("Fibroblast|Mesothelial", cc)]            <- "Fibroblasts"
  o[grepl("Mast", cc)]                              <- "Mast cells"
  o
}
lvcat <- c("T cells","Myeloid","NK cells","B cells","Plasma cells","Epithelial cells","Fibroblasts","Mast cells","Other")
mkset <- function(sub, lv) as.numeric(table(factor(sub, levels=lv)))/length(sub)
clong <- rbind(
  data.frame(group="all",   celltype=factor(lvcat,levels=lvcat), frac=mkset(mkcat(ct), lvcat)),
  data.frame(group="state", celltype=factor(lvcat,levels=lvcat), frac=mkset(mkcat(ct[SC$state]), lvcat)))
p_f <- ggplot(clong, aes(x=group, y=frac, fill=celltype)) +
  geom_bar(stat="identity", width=0.62) +
  scale_fill_manual(values=c(rev(brewer.pal(8,"Set2")),"#999999"), name=NULL) +
  scale_y_continuous(labels=function(x) sprintf("%.0f%%",x*100), breaks=seq(0,1,0.25)) +
  labs(title="Composition: state vs all cells", x=NULL, y="share of cells") +
  theme_classic(base_size=6.5, base_family="Arial") +
  theme(plot.title=element_text(size=8,face="bold",hjust=0), axis.title.y=element_text(size=6.5),
        axis.text=element_text(size=6.5), legend.position="right", legend.text=element_text(size=6),
        legend.key.size=unit(0.28,"cm"), plot.margin=margin(2,2,2,2)) +
  labs(tag="f")

## ---- panel g: double bars IGHG+ / FOS+ carrier vs rest (dodged) ----
pc <- read.csv(file.path(FIGS, "F12_plasma_caf_sub.csv")); pc$state <- as.logical(as.character(pc$state))
pgy <- pc |> group_by(lineage, sub) |> summarise(v=100*mean(state), n=n(), .groups="drop")
pl_ig  <- pgy$v[pgy$lineage=="Plasma cells" & pgy$sub==0]; pl_rest <- pgy$v[pgy$lineage=="Plasma cells" & pgy$sub==1]
fb_fos <- pgy$v[pgy$lineage=="Fibroblasts" & pgy$sub==1];  fb_rest <- pgy$v[pgy$lineage=="Fibroblasts" & pgy$sub==0]
pg <- data.frame(group=factor(c("IGHG+ plasma","IGHG+ plasma","FOS+ CAF","FOS+ CAF"),
                              levels=c("IGHG+ plasma","FOS+ CAF")),
                 type=c("carrier","rest","carrier","rest"), val=c(pl_ig,pl_rest,fb_fos,fb_rest))
p_g <- ggplot(pg, aes(x=group, y=val, fill=type)) +
  geom_col(position=position_dodge(0.62), width=0.5) +
  geom_text(aes(label=sprintf("%.1f",val)), position=position_dodge(0.62), vjust=-0.3, size=2.6, family="Arial") +
  scale_fill_manual(values=c(carrier=SIG, rest=GRA), name=NULL, labels=c("carrier","rest")) +
  scale_y_continuous(limits=c(0,7)) +
  labs(title="State fraction in carrier subsets", x=NULL, y="state fraction (%)") +
  theme_classic(base_size=6.5, base_family="Arial") +
  theme(plot.title=element_text(size=8,face="bold",hjust=0), axis.title.y=element_text(size=6.5),
        axis.text=element_text(size=6.5), legend.position="bottom", legend.text=element_text(size=6.5),
        legend.key.size=unit(0.22,"cm"), plot.margin=margin(2,2,2,2)) +
  labs(tag="g")

p_a <- panel_img(file.path(FIGS, "F2a_state_umap.png")) + labs(tag = "a")

## ---- panel h: carrier module structure -> two-arm diverging bars ----------------
arm <- data.frame(carrier=c("IGHG+ plasma","FOS+ CAF"), arm=c("M1 (immune)","M2/M3 (mito)"),
                  val=c(-0.106,-0.174), stringsAsFactors=FALSE)
p_h <- ggplot(arm, aes(x=carrier, colour=arm)) +
  geom_segment(aes(xend=carrier), y=-0.435, yend=-0.435, linewidth=0.35, colour="grey80") +
  geom_point(aes(y=val), size=2.8) +
  geom_segment(aes(y=0, yend=val), linewidth=1.2, lineend="butt") +
  geom_text(aes(y=val, label=sprintf("%.3f", val)), vjust=1.3, size=2.28, family="Arial") +
  scale_colour_manual(values=c("M1 (immune)"=SIG,"M2/M3 (mito)"=BLU), guide="none") +
  scale_y_continuous(limits=c(-0.50,0.06), breaks=c(-0.4,-0.2,0), labels=c("-0.4","-0.2","0")) +
  labs(x=NULL, y="module score", title="Carrier module structure") +
  theme(axis.text.x=element_text(size=6.5), axis.text.y=element_text(size=6.5),
        legend.position="none", plot.title=element_text(size=8)) +
  base_theme() + theme(plot.margin=margin(0,0,0,0,"mm")) + labs(tag="h")

## ---- panel i: module x set correlation matrix (RdBu CONTINUOUS) -------------------
icc <- data.frame(
  row=factor(c("M1\nimmune","M1\nimmune","M2\nmetab","M2\nmetab","M3\nmito-cyto","M3\nmito-cyto"),
             levels=c("M3\nmito-cyto","M2\nmetab","M1\nimmune")),
  col=factor(c("ICC45","random","ICC45","random","ICC45","random"), levels=c("ICC45","random")),
  r=c(0.635,0.004,0.700,0.004,0.700,0.004), stringsAsFactors=FALSE)
p_i <- ggplot(icc, aes(x=col, y=row, fill=r)) +
  geom_tile(colour="white", linewidth=0.5, height=1.0, width=0.95) +
  geom_text(aes(label=sprintf("%.2f", r)), size=2.28, colour="#111111", family="Arial", fontface="bold") +
  scale_fill_gradientn(colours=c("#2166AC","#F7F7F7","#B2182B"), values=c(0,0.0625,1),
                       limits=c(-0.05,0.75), name="Pearson r", breaks=c(0.0,0.70), labels=c("0.00","0.70")) +
  scale_x_discrete(position="top") +
  scale_y_discrete(expand=c(0.02,0.02)) +
  theme_void() +
  theme(axis.text.x.top=element_text(size=6.5, family="Arial", lineheight=0.9),
        axis.text.y=element_text(size=6.5, family="Arial", lineheight=0.9),
        legend.position="right", legend.text=element_text(size=6.5),
        legend.title=element_text(size=6.5), legend.key.width=unit(0.22,"cm"),
        legend.key.height=unit(0.62,"cm"), legend.margin=margin(0,0,0,0),
        plot.margin=margin(0,0,0,0,"mm"),
        plot.title=element_text(size=8, face="bold", hjust=0)) +
  labs(title="Module-set consistency (RdBu)", tag="i")

## =================================================================================
##  PART C -- panel j: pseudotime within carrier lineages, 2x2 paper tape
## =================================================================================
## tape sub-images are STRETCH-FILLED (crop ink-bbox then fill the cell) so the 2x2 tape
## has no internal blank bands; cell aspects ~= image aspects so stretch is minimal.
fill_panel <- function(path, bg="white", interpolate=TRUE) {
  img <- png::readPNG(path)
  bb <- img_bbox(img, 0.98); img <- img[bb$r0:bb$r1, bb$c0:bb$c1, , drop=FALSE]
  ggplot() +
    annotation_custom(grid::rasterGrob(img, interpolate=interpolate),
                      xmin=0, xmax=1, ymin=0, ymax=1) +
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    theme_void() +
    theme(plot.background=element_rect(fill=bg, colour=NA),
          panel.background=element_rect(fill=bg, colour=NA),
          plot.margin=margin(0,0,0,0))
}
p_jtl <- fill_panel(file.path(V9, "PT_fib_umap_traj.png"))
p_jtr <- fill_panel(file.path(V9, "PT_plasma_umap_traj.png"))
p_jbl <- fill_panel(file.path(V9, "PT_fib_state_along_c.png"))
p_jbr <- fill_panel(file.path(V9, "PT_plasma_state_along_c.png"))
p_grid <- (p_jtl | p_jtr) / (p_jbl | p_jbr) +
  plot_layout(widths=c(0.50,0.50), heights=c(0.56,0.44))
txt_row <- function(lab, size=2.28, face="plain", colour="black", bold=FALSE) {
  ggplot() + annotate("text", x=0, y=0.5, label=lab, hjust=0, size=size,
                      fontface=if (bold) "bold" else face, family="Arial", colour=colour) +
    scale_x_continuous(limits=c(-0.02,1), expand=c(0,0)) +
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) +
    theme_void() + theme(plot.margin=margin(0,0,0,0))
}
p_j <- (txt_row("pseudotime within carrier lineages (root: min-decouple cluster)", bold=TRUE) + labs(tag="j")) /
       p_grid /
       txt_row("epi lineage excluded (state 0.22%, too sparse)", size=2.12, colour="#555555") +
  plot_layout(heights=c(0.005, 0.991, 0.004))

## =================================================================================
##  PART D -- assemble v10.2 (vector panels; reading order a..j left->right per row).
##  Rows: a|b, c|d, e|f, g|h, i|j (j = full 2x2 pseudotime tape, widened to 3/4 width).
##  Row spacing 0.8mm; canvas 183mm; height content-compact.
## =================================================================================
design <- "AABB
CCDD
EEFF
GGHH
IJJJ"
fig <- p_a + p_b + p_c + p_d + p_e + p_f + p_g + p_h + p_i + p_j +
  plot_layout(design=design, heights=c(1.00,0.58,0.60,0.58,1.55),
              widths=c(1,1,1,1)) +
  tag_theme() +
  theme(plot.margin=margin(0.8,0.8,0.8,0.8,"mm"),
        plot.background=element_rect(fill="white", colour=NA))

save_pub_r(fig, "analysis/figs_final/_v9/v9_Fig2", width_mm=183, height_mm=215, dpi=600)
cat("Fig2 OK\n")
