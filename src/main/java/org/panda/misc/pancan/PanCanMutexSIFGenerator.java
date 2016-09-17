package org.panda.misc.pancan;

import org.panda.misc.MutexReader;
import org.panda.misc.MutexReader.Group;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

/**
 * For generating a in-same-group graph from Mutex results.
 */
public class PanCanMutexSIFGenerator
{
	public static void main(String[] args) throws IOException
	{
		String base = "/home/babur/Documents/PanCan/";
		int geneLimit = 500;
		prepare(base, geneLimit, base + "mutex-" + geneLimit);
	}

	private static void prepare(String dir, int geneLimit, String outSIFFile) throws IOException
	{
		Object[] o = PanCanResultLoader.readGroupsWithFlattenedControl(true, dir + "PanCan-results",
			dir + "PanCan-shuffled-?-results", f -> !f.getName().startsWith("random") && !f.getName().startsWith("PC"));

		Set<MutexReader.Group> groups = (Set<Group>) o[0];
		Map<String, Double> geneScores = MutexReader.convertGroupsToGeneBestScores(groups);

		double thr = geneScores.keySet().stream().sorted((g1, g2) -> geneScores.get(g1).compareTo(geneScores.get(g2)))
			.limit(geneLimit).mapToDouble(geneScores::get).max().getAsDouble();

		Map<String, Double> pairScores = MutexReader.convertGroupsToPairBestScores(
			groups.stream().filter(g -> g.score <= thr).collect(Collectors.toSet()));

		Set<String> genes = geneScores.keySet().stream().filter(g -> geneScores.get(g) <= thr).collect(Collectors.toSet());
		Map<String, Double> coverages = MutexToGraphAligner.readMutationCoverages(genes);

		ValToColor covColor = new ValToColor(new double[]{0, 0.2}, new Color[]{Color.WHITE, Color.RED});
		ValToColor sigColor = new ValToColor(new double[]{0, -Math.log(1E-3)}, new Color[]{Color.LIGHT_GRAY, Color.BLACK});

		BufferedWriter sifWriter = new BufferedWriter(new FileWriter(outSIFFile + ".sif"));
		BufferedWriter fmtWriter = new BufferedWriter(new FileWriter(outSIFFile + ".format"));

		fmtWriter.write("node\tall-nodes\tcolor\t255 255 255\n");
		fmtWriter.write("node\tall-nodes\tbordercolor\t0 0 0");

		pairScores.keySet().stream().forEach(k ->
		{
			String[] g = k.split(" ");
			String edge = g[0] + " in-same-group " + g[1];
			FileUtil.writeln(edge.replaceAll(" ", "\t"), sifWriter);
			FileUtil.lnwrite("edge\t" + edge + "\tcolor\t" + sigColor.getColorInString(-Math.log(pairScores.get(k))), fmtWriter);
		});

		geneScores.keySet().stream().filter(g -> geneScores.get(g) <= thr).forEach(g ->
		{
			FileUtil.lnwrite("node\t" + g + "\tcolor\t" + covColor.getColorInString(coverages.get(g)), fmtWriter);
			FileUtil.lnwrite("node\t" + g + "\ttooltip\t" + coverages.get(g), fmtWriter);
			FileUtil.lnwrite("node\t" + g + "\tbordercolor\t" + sigColor.getColorInString(-Math.log(geneScores.get(g))), fmtWriter);
		});

		sifWriter.close();
		fmtWriter.close();
	}
}
