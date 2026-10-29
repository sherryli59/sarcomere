// Actual Chromium/WebGL smoke check. Run on a compute node, never a login server.
const {chromium}=require(process.env.ENERGY_PLAYWRIGHT_MODULE);
const fs=require('node:fs');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader']});
 const page=await browser.newPage({viewport:{width:1550,height:1100},deviceScaleFactor:1});
 const errors=[],consoleErrors=[];
 page.on('pageerror',err=>errors.push(err.message));page.on('console',msg=>{if(msg.type()==='error')consoleErrors.push(msg.text());});
 await page.goto('file://'+process.argv[2],{timeout:120000,waitUntil:'load'});
 await page.waitForFunction(()=>document.body.dataset.ready==='true',{},{timeout:180000});
 await page.locator('#plot').screenshot({path:process.argv[3]+'/browser_whole.png'});
 await page.getByRole('button',{name:'+ Tip',exact:true}).click();
 await page.waitForFunction(()=>document.querySelector('#plot').layout.scene.xaxis.range[0]===580);
 await page.locator('#plot').screenshot({path:process.argv[3]+'/browser_tip.png'});
 await page.selectOption('#mode','1');
 await page.waitForFunction(()=>document.querySelector('#mode-note').textContent.includes('Integrated'));
 await page.check('#endpoint');await page.check('#slice');await page.uncheck('#shell-3');
 await page.waitForFunction(()=>document.querySelector('#plot').data.some(t=>t.type==='surface'));
 await page.fill('#px','700');await page.fill('#py','140');await page.locator('#pz').focus();
 await page.waitForFunction(()=>document.querySelector('#probe-readout').textContent.includes('Custom probe'));
 const probe=await page.evaluate(()=>window.energyViewer.evaluate([700,140,0],1));
 if(!(probe.energy>0))throw Error('Finite connector outside-midpoint example lost its contact');
 await page.locator('#plot').screenshot({path:process.argv[3]+'/browser_connector.png'});
 const result={passed:errors.length===0&&consoleErrors.length===0,page_errors:errors,console_errors:consoleErrors,
  outside_midpoint_connector:probe,webgl:await page.evaluate(()=>!!document.querySelector('#plot .gl-canvas')),
  checked:'Actual offline Chromium render; whole/tip views; finite mode; endpoint, slice and shell controls; custom probe'};
 fs.writeFileSync(process.argv[3]+'/browser_validation.json',JSON.stringify(result,null,2));
 if(!result.passed)throw Error(JSON.stringify(result));
 await browser.close();console.log('Actual browser checks passed');
})().catch(e=>{console.error(e);process.exitCode=1;});
