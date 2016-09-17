package org.panda.misc;

import org.panda.utility.ValToColor;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

/**
 * @author Ozgun Babur
 */
public class PrepareChiBESIFNodeColorFormatFile
{
	public static void prepare(Map<String, Double> valsMap, ValToColor vtc, String filename) throws IOException
	{
		BufferedWriter writer = new BufferedWriter(new FileWriter(filename));
		writer.write("node\tall-nodes\tcolor\t255 255 255");
		writer.write("\nnode\tall-nodes\tbordercolor\t0 0 0");
		for (String gene : valsMap.keySet())
		{
			writer.write("\nnode\t" + gene + "\tcolor\t" + vtc.getColorInString(valsMap.get(gene)));
			writer.write("\nnode\t" + gene + "\ttooltip\t" + valsMap.get(gene));
		}
		writer.close();
	}


}
