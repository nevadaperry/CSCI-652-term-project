// To run: npx ts-node this-file.ts

import * as fs from 'node:fs';
import * as readline from 'node:readline';
import { promisify } from 'node:util';

const BASES = ['A', 'C', 'G', 'T'] as const;
type Base = (typeof BASES)[number];
type Sequence = Base[];
const METABASES = [...BASES, '-', 'N'] as const;
type MetaBase = (typeof METABASES)[number];
type Matrix = number[][];
type Dendrogram = { [nodeName: string]: number };

main();

async function main(): Promise<void> {
	const maxVariants = Number(process.argv.at(-1)) || Number.MAX_SAFE_INTEGER;
	const path = 'provided/seqs';
	const seqNames = (
		(await promisify(fs.readdir)(path))
		// Exclude MERS and OC43 because they are weakly aligned outliers
		.filter(seqName => !['mers', 'oc43'].includes(seqName))
		.slice(0, maxVariants)
	);
	
	const legend = [...seqNames];
	const matrix = await initDistanceMatrix(seqNames);
	
	const dendrogram = await doUpgma(legend, matrix);
	console.log(`Final node heights: ${JSON.stringify(dendrogram, null, 2)}`);
}

async function initDistanceMatrix(seqNames: string[]): Promise<Matrix> {
	const matrix: Matrix =
		Array(seqNames.length)
		.fill(null)
		.map(() => Array(seqNames.length).fill(-1));
	for (let i = 0; i < seqNames.length; i++) {
		for (let j = i; j < seqNames.length; j++) {
			const score = await getScore(seqNames[i], seqNames[j]);
			matrix[i][j] = score;
			matrix[j][i] = score;
		}
	}
	return matrix;
}

async function getScore(seqName1: string, seqName2: string): Promise<number> {
	if (seqName1 == seqName2) return 0;
	let score: number;
	try {
		score = await readScoreFromMaf(seqName1, seqName2);
	} catch (e: any) {
		if (e.code !== 'ENOENT') throw e;
		try {
			score = await readScoreFromMaf(seqName2, seqName1);
		} catch (e: any) {
			if (e.code !== 'ENOENT') throw e;
			throw new Error(
				`File not found for ${seqName1} x ${seqName2}`
			);
		}
	}
	return score;
	async function readScoreFromMaf(
		seqName1: string,
		seqName2: string,
	): Promise<number> {
		const filePath = `provided/pw-more/${seqName1}.${seqName2}.sing.maf`;
		const readStream = fs.createReadStream(filePath);
		const readInterface = readline.createInterface(readStream);
		const readIterator = readInterface[Symbol.asyncIterator]();
		for (const skip of Array(11)) await readIterator.next();
		
		const { value }: { value: string } = await readIterator.next();
		if (!value.match('^a score=')) {
			throw new Error(`Unexpected line 12 in maf ${filePath}: ${value}`);
		}
		return Number(value.split('=')[1]);
	}
}

/**
 * Based on https://en.wikipedia.org/wiki/UPGMA#Working_example
 */
async function doUpgma(legend: string[], matrix: Matrix): Promise<Dendrogram> {
	const dendrogram: Dendrogram = {};
	while (true) {
		printStatus(legend, matrix);
		if (legend.length === 1) break;
		if (legend.length < 1) throw new Error();
		
		let lowestDist = Number.MAX_SAFE_INTEGER;
		let lowestI = -1, lowestJ = -1;
		for (let i = 0; i < matrix.length; i++) {
			for (let j = i + 1; j < matrix.length; j++) {
				if (matrix[i][j] < lowestDist) {
					lowestDist = matrix[i][j];
					lowestI = i;
					lowestJ = j;
				}
			}
		}
		const i = lowestI, j = lowestJ;
		const nodeName = `(${legend[j]}, ${legend[i]})`;
		const nodeHeight = Math.floor(matrix[i][j] / 2);
		dendrogram[nodeName] = nodeHeight;
		
		matrix.push(Array(legend.length + 1).fill(-1));
		for (let k = 0; k < legend.length; k++) {
			matrix[k].push(-1);
		}
		matrix[matrix.length - 1][matrix.length - 1] = 0;
		legend.push(nodeName);
		
		for (let k = 0; k < matrix.length - 1; k++) {
			if (k === i || k === j) continue;
			const combinedDist = (matrix[i][k] + matrix[j][k]) / 2;
			matrix[matrix.length - 1][k] = combinedDist;
			matrix[k][matrix.length - 1] = combinedDist;
		}
		removeRow(matrix, Math.max(i, j));
		removeRow(matrix, Math.min(i, j));
		removeCol(matrix, Math.max(i, j));
		removeCol(matrix, Math.min(i, j));
		legend.splice(Math.max(i, j), 1);
		legend.splice(Math.min(i, j), 1);
	}
	return dendrogram;
}

function removeRow(matrix: Matrix, rowNo: number) {
	matrix.splice(rowNo, 1);
}

function removeCol(matrix: Matrix, colNo: number) {
	for (let i = 0; i < matrix.length; i++) {
		matrix[i].splice(colNo, 1);
	}
}

function printStatus(legend: string[], matrix: Matrix) {
	console.log('Legend:');
	console.table(legend);
	console.log('Matrix:');
	console.table(matrix);
	console.log(``);
}

/* Unused by UPGMA
const filePaths = seqNames.map(seqName => `${path}/${seqName}`);
const seqs: Sequence[] = [];
for (const filePath of filePaths) {
	seqs.push(await getSeq(filePath));
	break;
}
async function getSeq(filePath: string): Promise<Sequence> {
	console.info(`Reading file ${filePath}`);
	const readStream = fs.createReadStream(filePath);
	const readInterface = readline.createInterface(readStream);
	const readIterator = readInterface[Symbol.asyncIterator]();
	
	// Skip first line because it's not part of the sequence
	readIterator.next();
	let seq = '';
	while (true) {
		const { done, value } = await readIterator.next();
		if (done) {
			console.info(`Reached end of file ${filePath}`);
			break;
		}
		seq += value;
	}
	return seq.trim() as unknown as Sequence;
}
*/
