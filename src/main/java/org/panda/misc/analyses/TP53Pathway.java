package org.panda.misc.analyses;

import org.panda.resource.network.SignedPCNoTransfac;
import org.panda.resource.signednetwork.SignedType;
import org.panda.utility.graph.DirectedGraph;
import org.pathwaycommons.sif.model.SIFGraph;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * @author Ozgun Babur
 */
public class TP53Pathway
{
	static SIFGraph load()
	{
		SIFGraph sifg = new SIFGraph();

		DirectedGraph graph = SignedPCNoTransfac.get().getGraph(SignedType.PHOSPHORYLATES);

		return null;
	}

	static void addRelations(DirectedGraph graph, SIFGraph sifg, String type)
	{
	}

	public static final Set<String> GENES = new HashSet<>(Arrays.asList("TP53 ATF3 ATM ATR AURKA BAK1 BRCA1 BRCA2 BRD7 CCNG1 CDK2 CDKN2A CHD8 CHEK1 CHEK2 CREBBP CSE1L CSNK1A1 CSNK1D CSNK1E DAXX DNMT3A DYRK2 E4F1 EP300 FBXO11 GPS2 GSK3B HIPK1 HIPK2 HIPK2 HUWE1 ING1 KAT2B KAT5 KAT8 KMT5A MAPK14 MAPK8 MAPK9 MDM2 MDM4 NEDD8 NOTCH1 PIN1 PLK1 PPM1D PPP1R13B PPP1R13L PPP2CA PRKCD PRMT5 PTEN PTPA RASSF1 RCHY1 RFWD2 RPL11 RPL23 RPL5 RPS6KA3 SETD2 SETD7 SKP2 SMARCA4 SMARCB1 SMYD2 STK11 TP53AIP1 TP53BP1 TP53BP2 TP53RK TRIM28 TTC5 UBE2D1 UBE3A USP7 VHL WT1 WWOX XPO1 YY1".split(" ")));
}
