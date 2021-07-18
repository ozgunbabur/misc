package org.panda.misc.analyses;

import org.panda.utility.ChartColorsList;
import org.panda.utility.FileUtil;
import org.panda.utility.TermCounter;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class COVID19
{
	public static final String DIR = "/Users/ozgun/Documents/Analyses/COVID19/";

	public static void main(String[] args) throws IOException
	{
		prepareFormatFile();
	}

	public static Map<String, Integer> getInteractionFreqs() throws IOException
	{
		TermCounter tc = new TermCounter();
		Files.lines(Paths.get(DIR + "SARS-Cov2-human-interacting-proteins.csv")).skip(2).map(l -> l.split("\t")[0].split(" ")[1]).forEach(tc::addTerm);
//		tc.print();
		return tc.getCountMap();
	}

	public static void prepareFormatFile() throws IOException
	{
		final String othersCol = "100 100 100";
		Map<String, Integer> freqs = getInteractionFreqs();

		List<String> vGenes = freqs.keySet().stream().sorted((o1, o2) -> freqs.get(o2).compareTo(freqs.get(o1))).collect(Collectors.toList());

		System.out.println("vGenes = " + vGenes);

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(DIR + "SARS-Cov2-interacting.format"));

		writer.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer.write("node\tall-nodes\tborderwidth\t2\n");

		int colSize = ChartColorsList.size();

		Files.lines(Paths.get(DIR + "SARS-Cov2-human-interacting-proteins.csv")).skip(2).map(l -> l.split("\t")).forEach(t ->
		{
			String hGene = t[2];
			String vGene = t[0].split(" ")[1];

			int index = vGenes.indexOf(vGene);
			if (index < colSize)
			{
				FileUtil.writeln("node\t" + hGene + "\tbordercolor\t" + ChartColorsList.getString(index), writer);
			}
			else
			{
				FileUtil.writeln("node\t" + hGene + "\tbordercolor\t" + othersCol, writer);
			}
			FileUtil.writeln("node\t" + hGene + "\ttooltip\t" + vGene, writer);
		});

		writer.close();

		BufferedWriter writer2 = Files.newBufferedWriter(Paths.get(DIR + "border-color-legend.sif"));
		vGenes.forEach(g -> FileUtil.writeln(g, writer2));
		writer2.close();
		BufferedWriter writer3 = Files.newBufferedWriter(Paths.get(DIR + "border-color-legend.format"));
		writer3.write("node\tall-nodes\tcolor\t255 255 255\n");
		writer3.write("node\tall-nodes\tborderwidth\t2\n");
		vGenes.forEach(g -> FileUtil.writeln("node\t" + g + "\tbordercolor\t" + (vGenes.indexOf(g) < colSize ? ChartColorsList.getString(vGenes.indexOf(g)) : othersCol), writer3));
		writer3.close();
	}
}
