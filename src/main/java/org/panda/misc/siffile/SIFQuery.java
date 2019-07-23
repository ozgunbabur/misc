package org.panda.misc.siffile;

import org.panda.resource.autismdatasets.SFARI;
import org.pathwaycommons.sif.model.*;
import org.pathwaycommons.sif.io.*;
import org.pathwaycommons.sif.query.*;
import org.pathwaycommons.sif.util.*;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static org.pathwaycommons.sif.util.EdgeAnnotationType.*;
import static org.pathwaycommons.sif.model.RelationTypeEnum.*;

public class SIFQuery
{
	public static void main(String[] args) throws IOException
	{
		findPathsBetween();
	}

	public static void findPathsBetween() throws IOException
	{
		// Tell the type and order of edge annotations in the SIF resource
		EdgeAnnotationType[] edgeAnnotTypes = new EdgeAnnotationType[]
			{
				DATA_SOURCE, PUBMED_IDS, PATHWAY_NAMES, MEDIATORS
			};

		// Initialize loader
		Loader loader = new Loader(edgeAnnotTypes);

		// Load a graph
		SIFGraph graph = loader.load(new FileInputStream("/home/ozgun/Data/PC/v10/PathwayCommons10.All.hgnc.txt"));

		//--- Perform a query on the graph

		// The query will traverse only two type of relations in this example
		EdgeSelector edgeSelector = new RelationTypeSelector(CONTROLS_STATE_CHANGE_OF, CONTROLS_EXPRESSION_OF);

		// Select a seed
		Set<String> seed = SFARI.get().getGenesWithMaxScore(6);
		System.out.println("seed = " + seed.size());

		// Run the query
		Set<Object> result = QueryExecutor.searchPathsBetween(graph, edgeSelector, seed, true, 1);

		//--- Report results

		edgeAnnotTypes = new EdgeAnnotationType[]
		{
			MEDIATORS
		};

		// Initialize the writer with the same edge annotation style
		Writer writer = new Writer(false, edgeAnnotTypes);

		FileOutputStream fos = new FileOutputStream("/home/ozgun/Documents/Temp/sif-query-result.sif");

		// Write results
		writer.write(result, fos);
	}


}
