# 扩展免疫 marker 库与 Stereo-seq RNA 可检测性审计：关键发现

## 1. Atlas 的免疫分类

Atlas 为 `BI_PF_ILD_atlas_v1`。共识别 9 个免疫 Level 1：B/Plasma; Basophil/Mast; Dendritic; Macrophage; Macrophage FABP4+; Monocyte; NK; T cell; pDC。
共识别 29 个免疫 Level 2：Activated NK XCL1+; B cell; B cell LRMP+; Basophil/Mast; CD8+ Naive T; Cytotoxic T; DC CXCL2+/CD14+/CD1C-; DC LAMP3+/CCR7+; DC1 CLEC9A+; DC2 C1QC+; DC2 CD1A+/CD1E+; DC2 CLEC10A+; Macrophage C1Q hi; Macrophage CHI3L1+/CD9 hi/; Macrophage FABP4+; Macrophage FABP4+/PDE4C+; Macrophage IL1B+; Macrophage LYVE1+; Macrophage RETN+/VCAN+; Monocytes CD14+; Monocytes CD14+/IL1B; Monocytes CD16+; NK; Plasma B; T CD4+ JUN+; T regulatory; helper T CD4+; helper T CD4+ CDR7+; pDC。Level 2 到 Level 1 的逐项核对见 `tables/02_ATLAS_CELLTYPE_HIERARCHY.tsv`。

## 2. 免疫 marker 数量

从免疫 Level 2 记录提取 1,450 条 marker 记录，涉及 692 个唯一基因。严格 lineage-core 为 27 个；进入扩展候选的非 core 基因为 209 个。

## 3. 谱系间重叠与特异性限制

有 181 个基因出现在两个或以上免疫 major lineages 的 atlas Top50 中；另有 40 个免疫 marker 同时出现在非免疫 Level 2 Top50 中。CD74/HLA-DRA、NKG7/GNLY、IGKC/IGHG 等均按 shared/functional 风险处理，不作为严格单一谱系证据。
本表只有各细胞类型的正向 Top marker，没有完整的非目标细胞表达矩阵。因此，“未在其他 Top50 出现”不等于真实细胞特异；本轮只能完成相对特异性审核，不能完成严格 sensitivity/specificity 验证。

## 4. Stereo-seq 中最容易检出的基因

MALAT1（21/21患者，84.523% Bin50）; TMSB4X（21/21患者，26.310% Bin50）; NTM（21/21患者，19.123% Bin50）; ACTB（21/21患者，18.663% Bin50）; CD74（21/21患者，17.630% Bin50）; DOCK4（21/21患者，17.455% Bin50）; PLXDC2（21/21患者，16.913% Bin50）; IGHGP（21/21患者，16.628% Bin50）; LPP（21/21患者，16.580% Bin50）; LSAMP（21/21患者，14.877% Bin50）; DPYD（21/21患者，14.747% Bin50）; TCF4（21/21患者，14.326% Bin50）; SMYD3（21/21患者，14.040% Bin50）; JMJD1C（21/21患者，13.835% Bin50）; SAT1（21/21患者，13.728% Bin50）。

## 5. 检出不足的经典 marker

CD79B（13/21患者，0.050% Bin50）; NKG7（16/21患者，0.078% Bin50）; IL1B（17/21患者，0.068% Bin50）; CD19（17/21患者，0.135% Bin50）; CD3D（18/21患者，0.072% Bin50）; KLRF1（18/21患者，0.113% Bin50）; CD3E（19/21患者，0.282% Bin50）; TRBC1（20/21患者，0.273% Bin50）; SPP1（20/21患者，0.736% Bin50）; S100A8（21/21患者，0.113% Bin50）; GNLY（21/21患者，0.166% Bin50）; S100A9（21/21患者，0.400% Bin50）。

## 6. 扩展 marker 是否增加覆盖

按八个 major lineages 分别计数并求和，原 small panels 在至少一位患者可检出 30 个 lineage-panel gene entries；extended panels 为 223 个，增加 193 个。该增加代表 detection coverage 增加，不代表免疫细胞识别准确率提高。每位患者的增量见 `tables/08_CORE_VS_EXTENDED_MARKER_COVERAGE.tsv`。

## 7. B / Plasma / T / NK / Myeloid 的 core 与 extended

- B core：BANK1;BLK;CD79B;MS4A1;TNFRSF13C；extended：ARHGAP24;CD22;ITSN2;KIAA0922;LINC00926;LRMP;MAP4K4;MGAT5;NCOA3;NLK;PRDM2;RBM6;SAMD12;TCL1A;USP34;VPREB3;ZCCHC7。
- Plasma core：无通过严格规则的基因；extended：DERL3;FAM46C;FCRL5;MZB1;PIM2;POU2AF1;SSR4;TXNDC5。
- T core：BCL11B;CD2;CD3D;CD3E;CD3G;ITK;TRAC;TRBC2；extended：ABCC1;AC016831.7;AC092580.4;ADAM19;ATXN1;CAMK4;CLEC2D;CTLA4;DUSP16;FAM129A;FOXO1;ICOS;IL2RA;PBX4;PDE3B;RNF19A;RP11-138A9.2;SPOCK2;THEMIS。
- NK core：无通过严格规则的基因；extended：CARD11;FGFBP2;KLRB1;SPON2;TXK;XCL1;XCL2。
- Macrophage core：FABP4;MARCO;MSR1;PPARG；extended：ABCA1;AC007192.4;ACP5;ALDH2;ANPEP;APOC1;APOE;ASAH1;ATP6V1F;BCAS4;CAPG;CCL20;CCSER1;CD163L1;CD52;CRIP1;CSTB;CTSB;CTSC;CTSL;CTSZ;CYP27A1;DOCK4;DOCK8;EMP3;EPB41L3;FABP5;FBP1;FMN1;FMNL2;FOLR2;FRMD4A;FTL;GLUL;GPNMB;GRN;H2AFY;HMOX1;LGALS1;LGALS3;LGMN;LILRB4;LIPA;LYVE1;MCEMP1;MS4A4A;MS4A6E;MT1G;NFKBIA;OLR1;PAPSS2;PLA2G7;PSAP;RNF130;RP11-20I23.1;RP11-295K3.1;SCHLAP1;SDC2;SH3BGRL3;SLC11A1;SLCO2B1;SNTB1;SNX10;TANC2;TREM1;TREM2;VSIG4;XXbac-BPG181M17.5。
- Monocyte core：CD300E;FCN1；extended：ADGRE2;CSF3R;FAM49A;IRAK3;LILRA5;LILRB2;LYN;MNDA;NAIP;PAG1;SLC25A37;SLC2A3;TNFRSF1B;WARS。
- DC core：CLEC10A;IL3RA；extended：ANKRD11;ANKRD33B;ARAP2;AXL;C1orf54;CCL22;CCR7;CD1A;CD1C;CD1E;CD207;CD86;CDYL;CERS6;CFLAR;CKLF;CLEC9A;CLIC2;CLN8;CPNE3;CPVL;CSF1R;CSF2RA;CST3;DAPP1;DUSP4;ENTPD1;ETV6;FAM160A1;FAM26F;FCER1A;FCGBP;FCGR2A;FCGR2B;FSCN1;GAB1;GNA15;GPR137B;GPR157;HDAC9;ID2;INPP4A;IRF4;IRF7;LGALS2;MCOLN2;NAAA;NAV1;NR4A3;PKIB;POLB;PTPN1;PTPRS;RAB31;RALA;REL;RFTN1;RGS10;RNASE6;RUBCN;S100B;SEPT6;SIPA1L3;SLC15A4;SLC22A23;SLC7A5;TBC1D8;WDFY4;YWHAH;ZBTB46。
- Mast/Basophil core：CPA3;HDC;IL1RL1;MS4A2;TPSAB1;TPSB2；extended：GATA2;LTC4S;MAOB;SLC18A2;SLC24A3;VWA5A。

## 8. 具备进一步探索条件的 Level 2 状态

Activated NK; Alveolar/resident macrophage; B; CD4 T; CD8 T; Cytotoxic T; Dendritic cell; Inflammatory macrophage; Mast/Basophil; Monocyte; NK; Plasma; Regulatory T; pDC
这些状态仅具备 RNA coverage 层面的候选条件，尚未证明在 50 µm neighborhood 中可被准确区分。

## 9. 可能无法稳定区分的状态

未观察到完全缺乏最低检测覆盖的 major state；仍需关注亚型间共享 marker。

## 10. 与 RCTD Reference18 的来源重叠

未确认。Excel 只给出 `BI_PF_ILD_atlas_v1` 名称和 marker 统计，没有 Reference18 构建来源、训练对象或样本 accession；因此不能证明两者独立，也不能证明存在来源重叠。这里只做标签对应，不把它当作独立验证。

## 11. 是否适合下一步 50 µm neighborhood scoring

可进入人工审核后的 exploratory scoring 准备，但必须同时保留三层：reference-based core、extended exploratory、Stereo-seq detectable subset，并把 shared/functional genes 单独呈现。不能因为扩展 panel 检出的基因更多，就声称细胞识别准确率提高；真正的准确率仍需独立参考、空间/蛋白证据或受控 benchmark。

## 本轮边界

未重新 RCTD、未修改 Reference18、未聚类、未改 K、未读取 niche enrichment 来优化 marker、未做 neighborhood scoring、DEG/GSEA、疾病检验或 CellChat。
