/* ============================================================
   文件名: global.js
   描述: 全站通用脚本 (包含 MathJax、底栏加载、Tab切换、复制功能等)
   适用页面: 首页、索引页、文章页
   ============================================================ */

/* ------------------------------------------------------------
   1. 全局配置 (必须在库加载前执行)
   ------------------------------------------------------------ */
window.MathJax = {
    tex: {
        inlineMath: [['$', '$'], ['\\(', '\\)']],
        displayMath: [['$$', '$$'], ['\\[', '\\]']]
    },
    options: {
        ignoreHtmlClass: 'tex2jax_ignore',
        processHtmlClass: 'tex2jax_process'
    }
};

/* ------------------------------------------------------------
   2. 全局工具函数 (供 HTML onclick 调用)
   ------------------------------------------------------------ */

   /* --- 全局变量：用于追踪当前的复制状态 --- */
   let copyFeedbackTimer = null;    // 记录当前的定时器 ID
   let currentActiveElement = null; // 记录当前正在显示“已复制”的元素
   let originalContentCache = "";   // 缓存该元素的原始 HTML 内容
   
   /**
    * 优化后的复制函数
    */
   function copyToClipboard(text, element, event) {
       // 1. 处理“点击另一个元素”的情况
       // 如果当前有元素正在显示反馈，且点击的是【另一个】元素
       if (currentActiveElement && currentActiveElement !== element) {
           resetCopyFeedback(); // 立即让上一个元素恢复原样
       }
   
       // 2. 处理“连续点击同一个元素”的情况
       // 如果点击的是同一个元素，我们需要清除之前的定时器，以便重新开始 1.5s 计时
       if (copyFeedbackTimer) {
           clearTimeout(copyFeedbackTimer);
       }
   
       // 3. 执行复制操作
       navigator.clipboard.writeText(text).then(() => {
           // 如果是第一次点击（或状态已复位），记录原始内容
           if (currentActiveElement !== element) {
               currentActiveElement = element;
               originalContentCache = element.innerHTML;
           }
   
           // 显示反馈样式
           element.innerHTML = '<i class="fas fa-check" style="color: #27ae60;"></i> 已复制';
           
           // 4. 设置 1.5秒后恢复
           // 无论是连续点击还是新点击，都会重新启动这 1.5秒的倒计时
           copyFeedbackTimer = setTimeout(() => {
               resetCopyFeedback();
           }, 1500);
   
       }).catch(err => {
           console.error('无法复制内容: ', err);
       });
   }
   
   /**
    * 内部辅助函数：将当前激活的元素恢复原状
    */
   function resetCopyFeedback() {
       if (currentActiveElement) {
           currentActiveElement.innerHTML = originalContentCache;
           currentActiveElement = null;
           originalContentCache = "";
       }
       if (copyFeedbackTimer) {
           clearTimeout(copyFeedbackTimer);
           copyFeedbackTimer = null;
       }
   }

/**
 * 集合切换函数 (用于索引页 Tab 切换)
 */
function showCollection(collectionId, btnElement) {
    // 1. 隐藏所有文章列表
    const lists = document.querySelectorAll('.article-list');
    lists.forEach(list => list.classList.remove('active'));
    
    // 2. 显示目标列表
    const target = document.getElementById(collectionId);
    if (target) target.classList.add('active');

    // 3. 按钮状态切换
    const btns = document.querySelectorAll('.collection-btn');
    btns.forEach(btn => btn.classList.remove('active'));
    if (btnElement) btnElement.classList.add('active');
}

/* ------------------------------------------------------------
   3. 内部辅助函数 (不直接暴露给 HTML)
   ------------------------------------------------------------ */

/**
 * 加载底栏
 */
function loadFooter() {
    // 优先读取页面中定义的 window.footerPath，如果没有定义，尝试读取 ./footer.html
    const path = window.footerPath || './footer.html';
    
    fetch(path)
        .then(response => {
            if (!response.ok) throw new Error(`Footer file not found at: ${path}`);
            return response.text();
        })
        .then(data => {
            const placeholder = document.getElementById('footer-placeholder');
            if (placeholder) placeholder.innerHTML = data;
        })
        .catch(err => console.error('底栏加载失败:', err));
}

/**
 * 强制预加载集合按钮的图标 (防止 Hover 闪烁)
 */
function preloadButtonImages() {
    const btns = document.querySelectorAll('.collection-btn');
    if (btns.length === 0) return;

    btns.forEach(btn => {
        const style = btn.getAttribute('style');
        if (!style) return;
        
        // 提取 CSS 变量中的 url路径
        const matches = style.match(/url\(['"]?([^'"]+)['"]?\)/g);
        if (matches) {
            matches.forEach(urlFunc => {
                const src = urlFunc.replace(/url\(['"]?|['"]?\)/g, '');
                new Image().src = src; // 触发浏览器缓存
            });
        }
    });
}

/**
 * 处理 URL 哈希跳转
 */
function handleHashNavigation() {
    const hash = window.location.hash.substring(1);
    if (hash) {
        const targetBtn = document.querySelector(`button[onclick*="'${hash}'"]`);
        if (targetBtn) {
            showCollection(hash, targetBtn);
        }
    }
}

/* ------------------------------------------------------------
   Prism.js 自动高亮加载器 (动态注入)
   ------------------------------------------------------------ */
function loadPrismHighlighter() {
    // 1. 检查当前页面是否有代码块，如果没有，直接退出，不浪费资源
    if (!document.querySelector('pre code')) return;

    // console.log("检测到代码块，正在加载 Prism 高亮引擎...");

    // 2. 定义资源路径 (使用稳定的 CDN)
    // 选用 "Tomorrow Night" 主题，完美匹配你的 #2c3e50 深色背景
    const cssUrl = 'https://cdnjs.cloudflare.com/ajax/libs/prism/1.29.0/themes/prism-tomorrow.min.css';
    const coreUrl = 'https://cdnjs.cloudflare.com/ajax/libs/prism/1.29.0/prism.min.js';
    const loaderUrl = 'https://cdnjs.cloudflare.com/ajax/libs/prism/1.29.0/plugins/autoloader/prism-autoloader.min.js';

    // 3. 动态加载 CSS
    const link = document.createElement('link');
    link.rel = 'stylesheet';
    link.href = cssUrl;
    document.head.appendChild(link);

    // 4. 动态加载 JS (核心库 -> 自动加载插件)
    // 辅助函数：加载脚本并返回 Promise
    const loadScript = (url) => {
        return new Promise((resolve, reject) => {
            const script = document.createElement('script');
            script.src = url;
            script.onload = resolve;
            script.onerror = reject;
            document.body.appendChild(script);
        });
    };

    // 按顺序加载：先加载核心，再加载自动识别插件
    loadScript(coreUrl)
        .then(() => loadScript(loaderUrl))
        .then(() => {
            // console.log("Prism 高亮引擎加载完成");
        })
        .catch(err => console.error("Prism 加载失败:", err));
}

/* ------------------------------------------------------------
   4. 页面初始化逻辑 (DOMContentLoaded)
   ------------------------------------------------------------ */
window.addEventListener('DOMContentLoaded', () => {
    
    // === A. 所有页面通用的逻辑 ===
    loadFooter(); // 加载底栏
    loadPrismHighlighter();

    // === B. 仅在“索引页”执行的逻辑 (通过 body class 判断) ===
    if (document.body.classList.contains('category-page')) {
        handleHashNavigation(); // 检查是否有 #rnaseq 这种哈希
    }

    // === C. 仅在“文章页”执行的逻辑 ===
    if (document.body.classList.contains('article-page')) {
        // 如果未来有文章页特有的逻辑(比如生成目录)，写在这里
    }
});