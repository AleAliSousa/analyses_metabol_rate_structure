// Run with the bundled artifact-tool runtime (see DATASET_REGION_COVERAGE.md).
import fs from 'node:fs/promises';
import path from 'node:path';
import { Workbook, SpreadsheetFile } from '@oai/artifact-tool';

const [input, outputDir] = process.argv.slice(2);
if (!input || !outputDir) throw new Error('Usage: node builder.mjs input.json output_directory');
const data = JSON.parse(await fs.readFile(input, 'utf8'));
await fs.mkdir(outputDir, { recursive: true });
const wb = Workbook.create();
const dark = '#273E55', muted = '#5D6B78';
const colors = {A:'#DAEBDD', S:'#DCE8F5', U:'#FFE4B7', P:'#E6E0F0', T:'#F7E5E8', '-':'#F3F4F5'};
const col = n => { let s=''; for(n++;n;n=Math.floor((n-1)/26))s=String.fromCharCode(65+(n-1)%26)+s; return s; };

function sheet(name, title, headers, rows, start, widths, height=38) {
  const s=wb.worksheets.add(name);
  s.showGridLines=false;
  const last=col(headers.length-1), end=start+rows.length;
  const all=s.getRange(`A1:${last}${end}`);
  all.format.font={name:'Arial',size:11,color:'#263442'};
  all.format.verticalAlignment='center';
  s.getRange('A2').values=[[title]];
  s.getRange('A2').format.font={size:16,bold:true,color:dark};
  s.getRange(`A3:${last}3`).format.borders={bottom:{style:'thin',color:'#BAC6D0'}};
  s.getRange(`A${start}:${last}${end}`).values=[headers,...rows];
  const table=s.tables.add(`A${start}:${last}${end}`,true,name.replaceAll(' ','')+'Table');
  table.showFilterButton=true;
  table.style='TableStyleLight1';
  s.getRange(`A${start+1}:${last}${end}`).format.fill='#FFFFFF';
  s.getRange(`A${start}:${last}${end}`).format.rowHeight=height;
  s.getRange(`A${start}:${last}${end}`).format.wrapText=true;
  s.getRange(`A${start}:${last}${start}`).format={fill:dark,font:{bold:true,color:'#FFFFFF'},horizontalAlignment:'center',rowHeight:42};
  widths.forEach((w,i)=>s.getRange(`${col(i)}1:${col(i)}${end}`).format.columnWidthPx=w);
  s.freezePanes.freezeRows(start);
  s.freezePanes.freezeColumns(name==='Coverage'?2:1);
  return s;
}

const rows=data.regions.map(({family,region})=>[family,region,...data.matrix.map(d=>d[region])]);
const cov=sheet('Coverage','Metabolic datasets: anatomical coverage',['Family','Anatomical territory',...data.datasets.map(d=>d.label)],rows,12,[155,300,145,145,145,155,155],28);
cov.tabColor=dark;
const context=[
  'Five audited papers. Shared names are candidate matches, not proof of identical anatomical boundaries.',
  '70 territories include nested and overlapping regions. They are not 70 independent anatomical observations.',
  'A  Mean available: named atlas ROI or cortical aggregate',
  'S  Mean available: selected spatial sample',
  'U  Mean available: label or tissue boundary needs clarification',
  'P  Published native means await transcription (Heiss 1984)',
  'T  Statistical tests only, no regional mean in the audited table',
  '–  No separately reported measurement in the audited tables',
];
context.forEach((v,i)=>{cov.getCell(i+3,0).values=[[v]];cov.getCell(i+3,0).format.font={color:i<2?muted:'#263442',size:11};});
cov.getRange('A4:G11').format.rowHeight=21;
cov.getRange(`C13:G${12+rows.length}`).format.horizontalAlignment='center';
for(const [code,fill] of Object.entries(colors)) {
  cov.getRange(`C13:G${12+rows.length}`).conditionalFormats.add('cellIs',{operator:'equal',formula:`"${code}"`,format:{fill,font:{color:code==='-'?'#A3AAB0':'#263442',bold:code!=='-'}}});
}
rows.forEach((r,i)=>r.slice(2).forEach((code,j)=>{
  cov.getCell(i+12,j+2).format={fill:colors[code],font:{color:code==='-'?'#A3AAB0':'#263442',bold:code!=='-'}};
}));
let previous='';
data.regions.forEach((r,i)=>{
  if(r.family!==previous)cov.getRange(`A${i+13}:G${i+13}`).format.borders={top:{style:'thin',color:'#9AAEBE'}};
  previous=r.family;
});

// One row per paper, with cohort/variant metadata beside the inventory.
const datasetHeaders=['Dataset','People (n)','Reported territories','Territories with available means','Population','Measure and units','Variants within cohort','Participant dependence','ROI method','Use in this project','Source','Limitations'];
const ds=sheet('Datasets','Source and cohort inventory',datasetHeaders,data.datasets.map(d=>[d.label,d.n,d.reported_territories,d.means_available_territories,d.population,`${d.measure}\n${d.units}`,d.variants,d.independence,d.method,d.role,d.source,d.caution]),5,[155,85,110,140,260,300,300,340,340,300,400,440],90);
ds.getRange('A4').values=[['People and regional coverage are different sample sizes. Publication counts do not establish cohort independence.']];
ds.getRange('A4').format.font={color:muted,size:11};
ds.getRange('B6:D10').setNumberFormat('0');

const evidenceHeaders=['Dataset','Family','Anatomical territory','Published label','Code','Mean availability','Anatomical scope','Anatomical / interpretation note','Source'];
const ev=sheet('Region evidence','Source-specific region evidence',evidenceHeaders,data.evidence.map(r=>[data.datasets.find(d=>d.dataset_id===r.dataset_id).label,r.family,r.region,r.source_label,r.coverage,r.mean_availability.replaceAll('_',' '),r.anatomical_scope.replaceAll('_',' '),r.notes,r.source]),5,[155,155,300,300,65,200,240,550,400],60);
ev.getRange('A4').values=[['Only reported territories appear here. The machine-readable evidence CSV also includes explicit absent cells.']];
ev.getRange('A4').format.font={color:muted,size:11};
ev.getRange(`E6:E${5+data.evidence.length}`).format.horizontalAlignment='center';
for(const [code,fill] of Object.entries(colors))ev.getRange(`E6:E${5+data.evidence.length}`).conditionalFormats.add('cellIs',{operator:'equal',formula:`"${code}"`,format:{fill}});
data.evidence.forEach((r,i)=>ev.getCell(i+5,4).format.fill=colors[r.coverage]);

wb.recalculate();
console.log((await wb.inspect({kind:'table',range:'Coverage!A12:G17',include:'values',tableMaxRows:6,tableMaxCols:7,maxChars:1800})).ndjson);
console.log((await wb.inspect({kind:'match',searchTerm:'#REF!|#DIV/0!|#VALUE!|#NAME\\?|#NUM!',options:{useRegex:true,maxResults:10},maxChars:1000})).ndjson);
for(const [name,range,suffix] of [
  ['Coverage','A1:G32','coverage_top'],['Coverage','A43:G70','coverage_subcortical'],
  ['Datasets','A1:F10','datasets_left'],['Datasets','G5:L10','datasets_right'],
  ['Region evidence','A1:E12','evidence_left'],['Region evidence','F5:I12','evidence_right'],
]) {
  const png=await wb.render({sheetName:name,range,scale:1,format:'png'});
  await fs.writeFile(path.join(outputDir,`${suffix}.png`),new Uint8Array(await png.arrayBuffer()));
}
const xlsx=await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(outputDir,'metabolic_dataset_region_coverage.xlsx'));
console.log('Exported workbook with 3 sheets.');
