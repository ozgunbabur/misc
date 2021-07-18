package org.panda.misc.causalpath;

import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Map;
import java.util.stream.Collectors;

public class CausalPathExternalAnnotationAdder
{
	static final Color maxDownColor = new Color(40, 80, 255);
	static final Color maxUpColor = new Color(255, 80, 40);

	public static void main(String[] args) throws IOException
	{
		printNewLines("/Users/ozgun/Downloads/p101_p3_diffs_May2021.txt", "x", "/Users/ozgun/Downloads/to-add-to-format-file.tsv");
//		temp();
	}

	private static void printNewLines(String valuesFile, String letter, String outFile) throws IOException
	{
		// read external numerical annotations
		Map<String, Double> valMap = FileUtil.lines(valuesFile).skip(1).map(l -> l.split("\t"))
			.collect(Collectors.toMap(t -> t[0], t -> Double.valueOf(t[1])));

		// find extremes
		double max = valMap.values().stream().max(Double::compareTo).get();
		double min = valMap.values().stream().min(Double::compareTo).get();


		ValToColor vtc = new ValToColor(new double[]{min, 0, max}, new Color[]{maxDownColor, Color.white, maxUpColor});

		BufferedWriter writer = Files.newBufferedWriter(Paths.get(outFile));

		valMap.keySet().stream().sorted().forEach(gene ->
		{
			double value = valMap.get(gene);
			FileUtil.writeln("node\t" + gene + "\trppasite\t" + gene + "-p3-score|" + letter + "|" + vtc.getColorInString(value) + "|100 100 100|" + value, writer);
		});

		writer.close();
	}

	public static void temp()
	{
		ValToColor vtc = new ValToColor(new double[]{-2, 0, 2}, new Color[]{maxDownColor, Color.white, maxUpColor});
		System.out.println("vtc.getColorInString(-1) = " + vtc.getColorInString(-1));
		System.out.println("vtc.getColorInString(1) = " + vtc.getColorInString(1));
	}
}
