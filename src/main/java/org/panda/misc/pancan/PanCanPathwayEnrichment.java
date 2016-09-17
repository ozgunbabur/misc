package org.panda.misc.pancan;

import org.biopax.paxtools.pattern.miner.SIFEnum;
import org.panda.misc.PathwayEnrichmentSIFGenerator;
import org.panda.resource.PCPathway;
import org.panda.utility.statistics.FDR;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.OptionalDouble;
import java.util.Set;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.toSet;

/**
 * @author Ozgun Babur
 */
public class PanCanPathwayEnrichment
{

	static void writePathwayEnrichment() throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		Set<String> genes = Files.lines(Paths.get(base + "pancan.txt")).skip(1).limit(300)
			.map(l -> l.split("\t")[0]).collect(Collectors.toSet());

		int minMember = 3;
		int maxMember = 300;
		PCPathway.get().writeEnrichmentResults(genes, minMember, maxMember, base + "pathway-enrichment.txt");
		Map<String, Double>[] pvals = PCPathway.get().getEnrichmentPvals(genes, null, minMember, maxMember);
		Map<String, Double> qvals = FDR.getQVals(pvals[0], pvals[1]);
		OptionalDouble thr = pvals[0].keySet().stream().filter(id -> qvals.get(id) < 0.1).mapToDouble(pvals[0]::get).max();

		Set<String> pathwayIDs = pvals[0].keySet().stream().filter(id -> pvals[0].get(id) <= thr.getAsDouble()).collect(toSet());

		PathwayEnrichmentSIFGenerator sg = new PathwayEnrichmentSIFGenerator();
		sg.setMolecules(genes);
		sg.showOnlySelectedMolecules(true);
		sg.setOwlFilename("/home/babur/Documents/PC/PathwayCommons.8.Detailed.BIOPAX.owl");
		sg.setBlacklistFilename("/home/babur/Documents/PC/blacklist.txt");
		sg.setPathwayIDs(pathwayIDs);
		String[][] groupTerms = {
			new String[]{"erbb", "egf"},
			new String[]{"c-met", "hgf"},
			new String[]{"shp2"},
			new String[]{"p53"},
			new String[]{"pdgf"},
			new String[]{"tgf_beta"},
			new String[]{"fgfr3"},
			new String[]{"insulin receptor"},
			new String[]{"estrogen receptor"},
			new String[]{"pi3 kinase"},
			new String[]{"jak-stat"},
			new String[]{"ras"},
//			new String[]{"miRNA"},
//			new String[]{"pyrimidine"},
//			new String[]{"Creatine"}
		};

		sg.setGroupTerms(groupTerms);
		sg.setTypes(SIFEnum.CONTROLS_STATE_CHANGE_OF, SIFEnum.CONTROLS_EXPRESSION_OF);

		sg.write(base + "pathway-enrichment");
	}
}
