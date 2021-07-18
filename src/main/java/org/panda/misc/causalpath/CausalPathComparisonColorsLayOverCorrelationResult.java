package org.panda.misc.causalpath;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

public class CausalPathComparisonColorsLayOverCorrelationResult
{
//	private static final String BASE = "/Users/ozgun/Documents/Analyses/CP-paper-runs/CPTAC-BRCA";
	private static final String BASE = "/Users/ozgun/Documents/Analyses/CPTAC-LUAD";
//	private static final String CORR_DIR = BASE + "/correlation-based-phospho";
	private static final String CORR_DIR = BASE + "/correlation-tumor";
//	private static final String COMP_DIR = BASE + "/Basal-vs-others";
	private static final String COMP_DIR = BASE + "/mutations/TP53";
	private static final String CORR_FILE = CORR_DIR + "/causative-apoptosis-neigh";

	private static final String CORR_SIF = CORR_FILE + ".sif";
	private static final String CORR_FORMAT = CORR_FILE + ".format";
	private static final String COMP_VALS = COMP_DIR + "/value-changes.txt";
	private static final String OUT_FORMAT = CORR_FILE + "-overlay-comp-mutTP53.format";

	public static void main(String[] args) throws IOException
	{
		run();
	}

	public static void run() throws IOException
	{
		Set<String> nodesInSIF = FileUtil.getTermsInTabDelimitedColumn(CORR_SIF, 0, 0);
		nodesInSIF.addAll(FileUtil.getTermsInTabDelimitedColumn(CORR_SIF, 2, 0));

//		Set<String> boxIDs = Files.lines(Paths.get(CORR_FORMAT)).map(l -> l.split("\t"))
//			.filter(t -> t[2].equals("rppasite"))
//			.map(t -> t[3].substring(0, t[3].indexOf("|"))).collect(Collectors.toSet());

		Map<String, Double> pMap = Files.lines(Paths.get(COMP_VALS)).map(l -> l.split("\t"))
			.filter(t -> t.length > 3 && !t[3].equals("NaN") && !t[3].equals("Q-value"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[2])));
		Map<String, Double> vMap = Files.lines(Paths.get(COMP_VALS)).map(l -> l.split("\t"))
			.filter(t -> t.length > 3 && !t[3].equals("NaN") && !t[3].equals("Q-value"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));

		ValToColor vtc = new ValToColor(new double[]{Math.log(0.0001), 0, -Math.log(0.0001)}, new Color[]{Color.BLUE, Color.WHITE, Color.RED});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(OUT_FORMAT));
		Files.lines(Paths.get(CORR_FORMAT)).forEach(l ->
		{
			String[] t = l.split("\t");

			if (t[0].equals("node") && t[2].equals("color") && vMap.containsKey(t[1]))
			{
				l = l.replace("255 255 255", getColor(t[1], vMap, pMap, vtc));
			}
			else if (t[2].equals("rppasite") && vMap.containsKey(t[3].substring(0, t[3].indexOf("|"))))
			{
				String id = t[3].substring(0, t[3].indexOf("|"));
				l = l.replace("255 255 255", getColor(id, vMap, pMap, vtc)) + ArrayUtil.getString(",", vMap.get(id), pMap.get(id), Math.signum(vMap.get(id)) * -Math.log(pMap.get(id)));
			}

			FileUtil.writeln(l, writer);
		});

//		boxIDs.stream().filter(vMap::containsKey)
//			.forEach(id -> FileUtil.writeln("node\t" + id + "\tcolor\t" +
//				getColor(id, vMap, pMap, vtc), writer));

		writer.close();
	}

	private static String getColor(String id, Map<String, Double> vMap, Map<String, Double> pMap, ValToColor vtc)
	{
		String white = "255 255 255";
		if (!vMap.containsKey(id)) return white;

		double v = vMap.get(id);
		double p = pMap.get(id);

		if (p >= 1) return white;

		return vtc.getColorInString(Math.signum(v) * -Math.log(p));
	}
}
